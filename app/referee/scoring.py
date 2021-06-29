import json
import os.path
import shutil

from django.db.models.fields.files import FieldFile

import ever_given.wrapper

from core import models


class AnswerPredictionPair:
    def __init__(self, answer, prediction):
        self.answer = answer
        self.prediction = prediction

    def __str__(self):
        return f"{self.answer} {self.prediction}"

    @staticmethod
    def args_from_keydict(keydict):
        args_dict = {}
        for key, pair in keydict.items():
            args_dict[f"{key}_answerkey"] = pair.answer
            args_dict[f"{key}_prediction"] = pair.prediction
        return args_dict


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

    scores_dict = utils.parse_output(scores_string)
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
