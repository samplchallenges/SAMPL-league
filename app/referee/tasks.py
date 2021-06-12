import json
import logging
import math
from collections import namedtuple

import dask
import dask.distributed as dd
import ever_given.wrapper

from core import models

logger = logging.getLogger(__name__)


AnswerPredictionPair = namedtuple("AnswerPredictionPair", ["answer", "prediction"])

SIMPLE = "SIMPLE"
MULTI = "MULTI_JSON"


def score_evaluation(container, evaluation, evaluation_score_types, output_handling):
    if len(evaluation_score_types) == 1:
        evaluation_score_handling = SIMPLE
        evaluation_score_type = list(evaluation_score_types.values())[0]
    else:
        evaluation_score_handling = MULTI

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

    command = None
    score_raw_args = {
        key: AnswerPredictionPair(answer_value, predictions[key])
        for key, answer_value in answer_keys.items()
    }

    if output_handling == SIMPLE:
        assert len(score_raw_args) == 1, "More than one output type for SIMPLE!"
        for ap_pair in score_raw_args.values():
            command = "{} {} ".format(ap_pair.answer, ap_pair.prediction)
    else:
        # TODO: this will generate JSON like
        # {"molwt": [78.4, 72], ..}
        # where the first list member is the answer key value
        # and the second is the prediction.
        # Perhaps a labeled dict structure would be better?
        command = json.dumps(score_raw_args)

    print(container.uri, command)
    scores_string = ever_given.wrapper.run_container(container.uri, command)
    if evaluation_score_handling == SIMPLE:
        score_value = float(scores_string)
        models.EvaluationScore.objects.create(
            evaluation=evaluation, score_type=evaluation_score_type, value=score_value
        )
    else:
        score_dict = json.loads(scores_string)
        for key, score_value in score_dict.items():
            models.EvaluationScore.objects.create(
                evaluation=evaluation,
                score_type=evaluation_score_types[key],
                value=score_value,
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
    if challenge.valuetype_set.filter(is_input_flag=False).count() == 1:
        output_handling = SIMPLE
    else:
        output_handling = MULTI

    evaluations = submission_run.evaluation_set.all()

    for evaluation in evaluations:
        score_evaluation(container, evaluation, evaluation_score_types, output_handling)

    run_scores_dicts = [
        {score.score_type.key: score.value for score in evaluation.scores.all()}
        for evaluation in evaluations
    ]

    command = json.dumps(run_scores_dicts)

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

    private_success.visualize(filename="delayed_graph.svg")
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


def _save_prediction(challenge, evaluation, output_type, raw_value):
    output_type_model = output_type.content_type.model_class()
    value_object = output_type_model.from_string(raw_value.strip())
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
    element = challenge.inputelement_set.get(pk=element_id, is_public=is_public)

    if challenge.valuetype_set.filter(is_input_flag=True).count() == 1:
        input_arg_handling = SIMPLE
    else:
        input_arg_handling = MULTI

    output_types = challenge.valuetype_set.filter(is_input_flag=False)
    if output_types.count() == 1:
        output_handling = SIMPLE
        output_type = output_types.first()
    else:
        output_handling = MULTI
        output_types_dict = {output_type.key: output_type
                             for output_type in output_types.all()}

    container = submission.container

    evaluation = models.Evaluation.objects.create(
        input_element=element, submission_run=submission_run
    )

    input_values_bykey = {
        input_value.value_type.key: input_value.value
        for input_value in element.inputvalue_set.all()
    }
    if input_arg_handling == SIMPLE:
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
        result = ever_given.wrapper.run_container(container.uri, command)
        if output_handling == SIMPLE:
            _save_prediction(challenge, evaluation, output_type, result)
        else:
            result_dict = json.loads(result)
            for key, value in result_dict.items():
                _save_prediction(challenge, evaluation, output_types[key], value)
        evaluation.status = models.Status.SUCCESS
    except:
        evaluation.status = models.Status.FAILURE
        raise
    finally:
        evaluation.save()
