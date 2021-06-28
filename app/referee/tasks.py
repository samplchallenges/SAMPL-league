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

from django.db.models.fields.files import FieldFile

import ever_given.wrapper

from core import models

logger = logging.getLogger(__name__)


SIMPLE = "SIMPLE"
MULTI = "MULTI"


class AnswerPredictionPair:
    def __init__(self, answer, prediction):
        self.answer = answer
        self.prediction = prediction

    def __str__(self):
        return f"{self.answer} {self.prediction}"

    @staticmethod
    def commandline_from_keydict(keydict):
        args_dict = {}
        for key, pair in keydict.items():
            args_dict[f"{key}_answerkey"] = pair.answer
            args_dict[f"{key}_prediction"] = pair.prediction
        return _prepare_commandline(args_dict)


def _prepare_commandline(args_dict):
    return " ".join(
        [f"--{key} {shlex.quote(value)}" for key, value in args_dict.items()]
    )


def _parse_output(raw_text):
    result = {}
    for line in raw_text.decode("utf-8").splitlines():
        lineparts = line.split(maxsplit=1)
        if len(lineparts) == 2:
            key, value = lineparts
            result[key] = value
        else:
            raise ValueError(f"Cannot parse output {line}, needs KEY VALUE format")
    return result


def _wrap_path(value, label):
    if isinstance(value, FieldFile):
        return os.path.join("/mnt", "inputs", label, os.path.basename(value.name))
    return value


def _wrap_prediction_path(value):
    return _wrap_path(value, "prediction")


def _wrap_answer_path(value):
    return _wrap_path(value, "answer")


def score_evaluation(container, evaluation, evaluation_score_types, scoringdir):

    predictions = {}
    for prediction in evaluation.prediction_set.all():
        predictions[prediction.value_type.key] = _wrap_prediction_path(prediction.value)
        if isinstance(prediction.value, FieldFile):
            shutil.copy(prediction.value.path, os.path.join("predictions", scoringdir))
    #TODO: make this faster by reusing outputdir for input

    input_element = evaluation.input_element


    answer_keys = {}
    for answer_key in input_element.answerkey_set.all():
        answer_keys[answer_key.value_type.key] = _wrap_answer_path(answer_key.value)
        if isinstance(answer_key.value, FieldFile):
            # Note: make this faster by only copying answer keys once per submission run
            shutil.copy(answer_key.value.path, os.path.join("answers", scoringdir))

    assert len(predictions) == len(answer_keys), "Error: number of predictions doesn't match answer keys, cannot score"
    score_raw_args = {
        key: AnswerPredictionPair(answer_value, predictions[key])
        for key, answer_value in answer_keys.items()
    }

    # This will generate a command line like
    # --molwt_answer 78.2 --molwt_prediction 72
    commandargs = AnswerPredictionPair.commandline_from_keydict(score_raw_args)
    command = f"score-evaluation {commandargs}"
    print(container.uri, command)
    scores_string = ever_given.wrapper.run_container(container.uri, command, scoringdir, None)

    scores_dict = _parse_output(scores_string)
    for key, score_value in scores_dict.items():
        models.EvaluationScore.objects.create(
            evaluation=evaluation,
            score_type=evaluation_score_types[key],
            value=float(score_value),
        )


def score_submission_run(container, submission_run, score_types):
    evaluation_score_types = score_types[models.ScoreType.Level.EVALUATION]
    submission_run_score_types = score_types[models.ScoreType.Level.SUBMISSION_RUN]

    if len(submission_run_score_types) == 1:
        submission_run_score_handling = SIMPLE
        submission_run_score_type = list(submission_run_score_types.values())[0]
    else:
        submission_run_score_handling = MULTI

    challenge = submission_run.submission.challenge

    evaluations = submission_run.evaluation_set.all()

    run_scores_dicts = [
        {score.score_type.key: score.value for score in evaluation.scores.all()}
        for evaluation in evaluations
    ]

    commandargs = json.dumps(run_scores_dicts)
    command = f"score-submissionrun {commandargs}"
    scores_string = ever_given.wrapper.run_container(container.uri, command)
    if submission_run_score_handling == SIMPLE:
        score_value = float(scores_string)
        models.SubmissionRunScore.objects.create(
            submission_run=submission_run,
            score_type=submission_run_score_type,
            value=score_value,
        )
    else:
        score_dict = json.loads(scores_string)
        for key, score_value in score_dict.items():
            models.SubmissionRunScore.objects.create(
                submission_run=submission_run,
                score_type=submission_run_score_types[key],
                value=score_value,
            )


def score_submission(submission_id, *run_ids):

    submission = models.Submission.objects.get(pk=submission_id)
    challenge = submission.challenge
    container = models.ScoreMaker.objects.get(challenge=challenge).container

    score_types = {
        models.ScoreType.Level.EVALUATION: {},
        models.ScoreType.Level.SUBMISSION_RUN: {},
    }

    for score_type in challenge.scoretype_set.all():
        score_types[score_type.level][score_type.key] = score_type

    for run_id in run_ids:

        submission_run = submission.submissionrun_set.get(pk=run_id)

        score_submission_run(container, submission_run, score_types)


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
    output_type_model = output_type.content_type.model_class()
    value_object = output_type_model.from_string(raw_value.strip(), output_dir=output_dir)
    value_object.save()
    models.Prediction.objects.create(
        challenge=challenge,
        value_type=output_type,
        evaluation=evaluation,
        value_object=value_object,
    )


@dask.delayed(pure=False)
def run_element(submission_id, element_id, submission_run_id, is_public):
    submission = models.Submission.objects.get(pk=submission_id)
    challenge = submission.challenge
    submission_run = submission.submissionrun_set.get(pk=submission_run_id)
    evaluation_score_types = {score_type.key: score_type for score_type in  challenge.scoretype_set.filter(level=models.ScoreType.Level.EVALUATION)}

    element = challenge.inputelement_set.get(pk=element_id, is_public=is_public)

    if challenge.valuetype_set.filter(is_input_flag=True).count() == 1:
        input_arg_handling = SIMPLE
    else:
        input_arg_handling = MULTI

    blob_content_type = models.ContentType.objects.get_for_model(models.BlobValue)
    file_content_type = models.ContentType.objects.get_for_model(models.FileValue)

    output_types = challenge.valuetype_set.filter(is_input_flag=False)
    output_types_dict = {
        output_type.key: output_type for output_type in output_types.all()
    }
    has_output_files = output_types.filter(content_type__in=(blob_content_type, file_content_type)).exists()

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
                value_type__content_type__in=(blob_content_type, file_content_type)
            )
        }

        for input_value in element.inputvalue_set.filter(
            value_type__content_type=blob_content_type
        ):
            filename = input_value.value_type.description  # TODO: hacky overlay
            if not filename:
                raise ValueError("description must be set to filename on blob types")
            input_values_bykey[input_value.value_type.key] = os.path.join(
                "/mnt", "inputs", filename
            )
            with open(os.path.join(inputdir, filename), "wb") as fp:
                fp.write(input_value.value)
        for input_value in element.inputvalue_set.filter(
            value_type__content_type=file_content_type
        ):
            filename = input_value.value.name
            path = input_value.value.path
            input_values_bykey[input_value.value_type.key] = os.path.join(
                "/mnt", "inputs", os.path.basename(filename)
            )
            shutil.copy(path, inputdir)
        command = _prepare_commandline(input_values_bykey)
        print(command)

        try:
            result = ever_given.wrapper.run_container(
                container.uri, command, inputdir=inputdir, outputdir=outputdir
            )
            result_dict = _parse_output(result)
            for key, value in result_dict.items():
                output_type = output_types_dict.get(key)
                if output_type:
                    _save_prediction(challenge, evaluation, output_type, value, outputdir)
                else:
                    print(f"Ignoring key {key} with value {value}")
            score_evaluation(challenge.scoremaker.container, evaluation, evaluation_score_types, scoringdir)
            evaluation.status = models.Status.SUCCESS
        except:
            evaluation.status = models.Status.FAILURE
            raise
        finally:
            evaluation.save()
