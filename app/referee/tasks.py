import json
import logging
import math

import dask
import dask.distributed as dd

import ever_given.wrapper

from core import models


logger = logging.getLogger(__name__)


def score_submission(submission_id, *run_ids):

    submission = models.Submission.objects.get(pk=submission_id)
    challenge = submission.challenge

    output_types = challenge.valuetype_set.filter(is_input_flag=False)
    # container = models.ScoreMaker.objects.get(challenge=challenge).container
    # container_uri = f"{container.registry}/{container.label}:{container.tag}"

    for run_id in run_ids:

        submission_run = submission.submissionrun_set.get(pk=run_id)
        evaluations = submission_run.evaluation_set.all()
        score_type, _ = models.ScoreType.objects.get_or_create(
            challenge=challenge, key="diff"
        )
        for evaluation in evaluations:

            predictions = {
                prediction.value_type.key: prediction.value
                for prediction in evaluation.prediction_set.all()
            }
            input_element = evaluation.input_element

            answer_keys = {
                answer_key.value_type.key: answer_key.value
                for answer_key in input_element.answerkey_set.all()
            }

            # if logger.isDebugEnabled():
            # smiles_type = challenge.valuetype_set.get(key="SMILES")
            # input_value = input_element.inputvalue_set.first()
            # get(input_type=smiles_type)
            #    logger.debug(
            #        "Prediction: %s %s answer: %s",
            #        prediction.value,
            #        input_value.value,
            #        answer_key.value,
            #    )
            # command = "{} {}".format(answer_key.value, prediction.value)
            # TODO: actually run the submission container
            # result = ever_given.wrapper.run_container(container_uri, command)
            score_sq = 0
            for key, answer_value in answer_keys.items():
                prediction_value = predictions.get(key)
                if not prediction_value:
                    raise ValueError(
                        "No prediction for key {} on evaluation {}".format(
                            key, evaluation
                        )
                    )
                score_sq += pow(answer_value - prediction_value, 2)
            score = math.sqrt(score_sq / len(answer_keys))
            models.EvaluationScore.objects.create(
                score_type=score_type, evaluation=evaluation, value=score
            )


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


@dask.delayed(pure=False)
def run_element(submission_id, element_id, submission_run_id, is_public):
    submission = models.Submission.objects.get(pk=submission_id)
    challenge = submission.challenge
    submission_run = submission.submissionrun_set.get(pk=submission_run_id)
    element = challenge.inputelement_set.get(pk=element_id, is_public=is_public)

    output_type = models.ValueType.objects.get(challenge=challenge, is_input_flag=False)
    container = submission.container

    evaluation = models.Evaluation.objects.create(
        input_element=element, submission_run=submission_run
    )
    container_uri = f"{container.registry}/{container.label}:{container.tag}"
    input_values_bykey = {
        input_value.value_type.key: input_value.value
        for input_value in element.inputvalue_set.all()
    }
    if len(input_values_bykey) == 1:
        input_arg = list(input_values_bykey.values())[0]
    else:
        input_arg = json.dumps(input_values_bykey)
    try:
        command = submission.challenge.execution_options_json["command"]
    except KeyError:  # if no execution options nothing to prepend
        command = ""
    command += input_arg
    print(command)

    try:
        result = ever_given.wrapper.run_container(container_uri, command)
        output_type_model = output_type.content_type.model_class()

        if output_type_model == models.FloatValue:
            result = float(result.strip())
            result_obj = models.FloatValue.objects.create(value=result)
        else:
            raise Exception("must only float value for output now")
        print(result)
        # maybe we return the prediction so we can score it?
        prediction = models.Prediction.objects.create(
            challenge=challenge,
            value_type=output_type,
            evaluation=evaluation,
            value_object=result_obj,
        )
        evaluation.status = models.Status.SUCCESS
        return prediction.pk
    except:
        evaluation.status = models.Status.FAILURE
        raise
    finally:
        evaluation.save()
