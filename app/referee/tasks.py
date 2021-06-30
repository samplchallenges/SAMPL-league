import json
import logging
import math
import os
import shlex
import shutil
import tempfile
from collections import namedtuple

import dask
import dask.distributed as dd

import ever_given.wrapper

from core import models

from .scoring import score_evaluation, score_submission_run


logger = logging.getLogger(__name__)


SIMPLE = "SIMPLE"
MULTI = "MULTI"


def run_and_score_submission(client, submission):
    """
    Runs public and private, plus scoring
    """
    challenge = submission.challenge
    public_element_ids = challenge.inputelement_set.filter(is_public=True).values_list(
        "id", flat=True
    )
    private_element_ids = challenge.inputelement_set.filter(
        is_public=False
    ).values_list("id", flat=True)

    public_run_id, public_prediction_ids = run_submission(
        submission.pk, public_element_ids, True, is_public=True
    )
    public_success = check_and_score(public_run_id, public_prediction_ids)
    private_run_id, private_prediction_ids = run_submission(
        submission.pk, private_element_ids, public_success, is_public=False
    )

    private_success = check_and_score(private_run_id, private_prediction_ids)

    # private_success.visualize(filename="delayed_graph.svg")
    future = client.submit(private_success.compute)
    print("Future key:", future.key)

    dd.fire_and_forget(future)
    return future


@dask.delayed(pure=False)
def check_and_score(submission_run_id, prediction_ids):
    submission_run = models.SubmissionRun.objects.get(pk=submission_run_id)
    submission_run.status = models.Status.SUCCESS
    challenge = submission_run.submission.challenge
    submission_run.save()

    print(
        "Running check_and_score",
        submission_run_id,
        "public?",
        submission_run.is_public,
    )
    score_submission(submission_run.submission.pk, submission_run_id)
    return True


@dask.delayed(pure=False)
def chain(upstream_success, delayed_func, *args, **kwargs):
    # Can ignore upstream_success
    return delayed_func(*args, **kwargs)


@dask.delayed(pure=False)
def create_submission_run(submission_id, conditional, is_public=True):
    # conditional will be a dask delayed; if it's false, the run_element will no-op
    if not conditional:
        return
    submission = models.Submission.objects.get(pk=submission_id)
    container = submission.container
    if not container.digest:
        container.digest = "DEADBEEF"
        container.save()
    submission_run = models.SubmissionRun.objects.create(
        submission=submission,
        digest=container.digest,
        is_public=is_public,
        status=models.Status.PENDING,
    )
    # TODO: need to store future key?
    # submission run pair is place to store?
    # submission_run.digest = future.key
    submission_run.save()
    return submission_run.id


def run_submission(submission_id, element_ids, conditional, is_public=True):

    submission_run_id = create_submission_run(
        submission_id, conditional, is_public=is_public
    )
    delayeds = dask.delayed(
        [
            run_element(
                submission_id,
                element_id,
                submission_run_id,
                is_public=is_public,
            )
            for element_id in element_ids
        ],
        nout=len(element_ids),
    )
    return (submission_run_id, delayeds)


def _save_prediction(challenge, evaluation, output_type, raw_value, output_dir):
    prediction = models.Prediction(
        challenge=challenge,
        value_type=output_type,
        evaluation=evaluation,
    )

    output_type_model = output_type.content_type.model_class()
    value_object = output_type_model.from_string(
        raw_value.strip(), prediction=prediction, output_dir=output_dir
    )
    value_object.save()


@dask.delayed(pure=False)
def run_element(submission_id, element_id, submission_run_id, is_public):
    submission = models.Submission.objects.get(pk=submission_id)
    challenge = submission.challenge
    submission_run = submission.submissionrun_set.get(pk=submission_run_id)
    evaluation_score_types = {
        score_type.key: score_type
        for score_type in challenge.scoretype_set.filter(
            level=models.ScoreType.Level.EVALUATION
        )
    }

    element = challenge.inputelement_set.get(pk=element_id, is_public=is_public)

    if challenge.valuetype_set.filter(is_input_flag=True).count() == 1:
        input_arg_handling = SIMPLE
    else:
        input_arg_handling = MULTI

    file_content_type = models.ContentType.objects.get_for_model(models.FileValue)

    output_types = challenge.valuetype_set.filter(is_input_flag=False)
    output_types_dict = {
        output_type.key: output_type for output_type in output_types.all()
    }
    has_output_files = output_types.filter(
        content_type__in=(blob_content_type, file_content_type)
    ).exists()

    container = submission.container

    evaluation = models.Evaluation.objects.create(
        input_element=element, submission_run=submission_run
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        inputdir = os.path.join(str(tmpdir), "input")
        os.mkdir(inputdir)
        if has_output_files:
            outputdir = os.path.join(str(tmpdir), "output")
            os.mkdir(outputdir)
            scoringdir = os.path.join(str(tmpdir), "scoring")
            os.mkdir(scoringdir)
            os.mkdir(os.path.join(scoringdir, "predictions"))
            os.mkdir(os.path.join(scoringdir, "answers"))
        else:
            outputdir = None
            scoringdir = None
        input_values_bykey = {
            input_value.value_type.key: input_value.value
            for input_value in element.inputvalue_set.exclude(
                value_type__content_type=file_content_type
            )
        }

        for input_value in element.inputvalue_set.filter(
            value_type__content_type=file_content_type
        ):
            filename = input_value.value.name
            path = input_value.value.path
            input_values_bykey[input_value.value_type.key] = os.path.join(
                "/mnt", "inputs", os.path.basename(filename)
            )
            shutil.copy(path, inputdir)
        command = utils.prepare_commandline(input_values_bykey)
        print(command)

        try:
            result = ever_given.wrapper.run_container(
                container.uri, command, inputdir=inputdir, outputdir=outputdir
            )
            result_dict = utils.parse_output(result)
            for key, value in result_dict.items():
                output_type = output_types_dict.get(key)
                if output_type:
                    _save_prediction(
                        challenge, evaluation, output_type, value, outputdir
                    )
                else:
                    print(f"Ignoring key {key} with value {value}")
            score_evaluation(
                challenge.scoremaker.container,
                evaluation,
                evaluation_score_types,
                scoringdir,
            )
            evaluation.status = models.Status.SUCCESS
        except:
            evaluation.status = models.Status.FAILURE
            raise
        finally:
            evaluation.save()
