import json
import logging
import tempfile

import ever_given.wrapper

from core import models

logger = logging.getLogger(__name__)


class MissingKeyException(Exception):
    pass


class ScoringFailureException(Exception):
    pass


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


def _build_kwargs(evaluation):
    predictions_by_key, file_predictions_by_key = models.Prediction.dicts_by_key(
        evaluation.prediction_set.all()
    )
    input_element = evaluation.input_element

    answer_keys, file_answer_keys = models.AnswerKey.dicts_by_key(
        input_element.answerkey_set.all()
    )

    assert len(predictions_by_key) == len(answer_keys), (
        f"Error: number of predictions ({len(predictions_by_key)}) "
        f"doesn't match answer keys ({len(answer_keys)}) , cannot score\n"
    )

    assert len(file_predictions_by_key) == len(file_answer_keys), (
        f"Error: number of file predictions ({len(file_predictions_by_key)}) "
        f"doesn't match answer keys ({len(file_answer_keys)}) , cannot score\n"
    )

    score_args = {
        key: AnswerPredictionPair(answer_value, predictions_by_key[key])
        for key, answer_value in answer_keys.items()
    }
    score_file_args = {
        key: AnswerPredictionPair(answer_value, file_predictions_by_key[key])
        for key, answer_value in file_answer_keys.items()
    }

    kwargs = AnswerPredictionPair.args_from_keydict(score_args)
    file_kwargs = AnswerPredictionPair.args_from_keydict(score_file_args)
    return kwargs, file_kwargs


def score_evaluation(container, evaluation, evaluation_score_types):
    kwargs, file_kwargs = _build_kwargs(evaluation)
    command = "score-evaluation"
    evaluation.append(stdout=f"Scoring with {container.uri} {command}\n")
    kwargs.update(container.custom_args())
    file_kwargs.update(container.custom_file_args())

    try:
        for key, score_value in ever_given.wrapper.run(
            container.uri,
            command,
            file_kwargs=file_kwargs,
            kwargs=kwargs,
            log_handler=models.Evaluation.LogHandler(evaluation),
        ):
            if key in evaluation_score_types:
                models.EvaluationScore.objects.create(
                    evaluation=evaluation,
                    score_type=evaluation_score_types[key],
                    value=float(score_value),
                )
    except Exception as exc:  # pylint: disable=broad-except
        raise ScoringFailureException from exc


def _score_submission_run(container, submission_run, score_types):
    submission_run_score_types = score_types[models.ScoreType.Level.SUBMISSION_RUN]

    evaluations = submission_run.evaluation_set.all()

    run_scores_dicts = [
        {score.score_type.key: score.value for score in evaluation.scores.all()}
        for evaluation in evaluations
    ]

    with tempfile.NamedTemporaryFile(suffix=".json", mode="w") as fp:
        json.dump(run_scores_dicts, fp)
        fp.flush()
        kwargs = {}
        file_kwargs = {"scores": fp.name}
        kwargs.update(container.custom_args())
        file_kwargs.update(container.custom_file_args())
        command = "score-submissionrun"
        matched_keys = set()
        for key, value in ever_given.wrapper.run(
            container.uri, command, file_kwargs=file_kwargs, kwargs=kwargs
        ):
            if key in submission_run_score_types:
                models.SubmissionRunScore.objects.create(
                    submission_run=submission_run,
                    score_type=submission_run_score_types[key],
                    value=value,
                )
                matched_keys.add(key)

        if matched_keys != set(submission_run_score_types.keys()):
            error_message = (
                f"Scoring submission run failed. Found keys of "
                f"{matched_keys} out of required {submission_run_score_types}"
            )
            raise MissingKeyException(error_message)


def score_submission_run(submission_run):
    submission = submission_run.submission
    challenge = submission.challenge
    container = models.ScoreMaker.objects.get(challenge=challenge).container

    try:
        _score_submission_run(container, submission_run, challenge.score_types)
    except Exception as exc:
        submission_run.append(stderr=str(exc))
        submission_run.status = models.Status.FAILURE
        raise ScoringFailureException from exc
    finally:
        submission_run.save(update_fields=["status"])
    submission_run.append(stdout="Scoring complete")
