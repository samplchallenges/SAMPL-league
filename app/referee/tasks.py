import logging
import tempfile
from pathlib import Path

import dask
import dask.distributed as dd
import ever_given.wrapper
from django.conf import settings
from ever_given.log_processing import CancelledException

from core import models

from . import scoring

logger = logging.getLogger(__name__)


def run_and_score_submission(client, submission):
    """
    Runs public and private, plus scoring
    """
    delayed_conditional = dask.delayed(True)
    for is_public in (True, False):
        delayed_conditional = _trigger_submission_run(
            submission, delayed_conditional, is_public=is_public
        )

    if settings.VISUALIZE_DASK_GRAPH:
        delayed_conditional.visualize(filename="task_graph.svg")

    future = client.submit(delayed_conditional.compute)  # pylint:disable=no-member
    logger.info("Future key: %s", future.key)

    dd.fire_and_forget(future)
    return future


def _trigger_submission_run(submission, delayed_conditional, *, is_public):
    submission_run = submission.create_run(is_public=is_public)
    evaluation_statuses = _run_evaluations(submission_run, delayed_conditional)
    return check_and_score(submission_run.id, delayed_conditional, evaluation_statuses)


@dask.delayed(pure=False)  # pylint:disable=no-value-for-parameter
def check_and_score(submission_run_id, delayed_conditional, evaluation_statuses):
    submission_run = models.SubmissionRun.objects.get(pk=submission_run_id)
    uniq_statuses = set(evaluation_statuses)
    if not delayed_conditional:
        status = models.Status.CANCELLED
    elif {models.Status.PENDING, models.Status.RUNNING} & uniq_statuses:
        submission_run.append(
            stderr=f"Evaluations should have all completed, but have statuses {evaluation_statuses}!"
        )
        status = models.Status.FAILURE
    elif {models.Status.CANCELLED} == uniq_statuses:
        status = models.Status.CANCELLED
    elif {models.Status.FAILURE, models.Status.CANCELLED} & uniq_statuses:
        status = models.Status.FAILURE
    else:
        status = models.Status.SUCCESS

    submission_run.status = status
    if status != models.Status.SUCCESS:
        submission_run.append(stderr=f"Submission run failed {status}")
        submission_run.save(update_fields=["status"])
        return False
    submission_run.append(stdout="Running check_and_score")
    scoring.score_submission_run(submission_run)

    return True


def _run_evaluations(submission_run, conditional):

    return [
        run_evaluation(
            submission_run.submission.id,
            evaluation.id,
            submission_run.id,
            conditional=conditional,
        )
        for evaluation in submission_run.evaluation_set.all()
    ]


@dask.delayed(pure=False)  # pylint:disable=no-value-for-parameter
def run_evaluation(submission_id, evaluation_id, submission_run_id, conditional):
    submission = models.Submission.objects.get(pk=submission_id)
    container = submission.container
    challenge = submission.challenge
    submission_run = submission.submissionrun_set.get(pk=submission_run_id)

    if not conditional or submission_run.check_cancel_requested():
        models.Evaluation.objects.filter(pk=evaluation_id).update(
            status=models.Status.CANCELLED
        )
        return models.Status.CANCELLED

    evaluation_score_types = challenge.score_types[models.ScoreType.Level.EVALUATION]

    evaluation = submission_run.evaluation_set.get(pk=evaluation_id)
    element = evaluation.input_element

    output_file_keys = challenge.output_file_keys()

    kwargs, file_kwargs = element.all_values()
    evaluation.mark_started(kwargs, file_kwargs)
    kwargs.update(container.custom_args())
    file_kwargs.update(container.custom_file_args())
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            dirpath = Path(str(tmpdir))
            output_dir = None
            if output_file_keys:
                output_dir = dirpath / "output"
                output_dir.mkdir()
            parsed_results = ever_given.wrapper.run(
                container.uri,
                kwargs=kwargs,
                file_kwargs=file_kwargs,
                output_dir=output_dir,
                output_file_keys=output_file_keys,
                log_handler=models.Evaluation.LogHandler(evaluation),
                cancel_requested_func=submission_run.check_cancel_requested,
                aws_login_func=settings.AWS_LOGIN_FUNCTION
                if settings.LOGIN_TO_AWS
                else None,
            )

            for key, value in parsed_results:
                output_type = challenge.output_type(key)
                if output_type:
                    prediction = models.Prediction.load_output(
                        challenge, evaluation, output_type, value
                    )
                    prediction.save()
                else:
                    evaluation.append(stderr=f"Ignoring key {key} with value {value}\n")

        scoring.score_evaluation(
            challenge.scoremaker.container,
            evaluation,
            evaluation_score_types,
        )
        evaluation.status = models.Status.SUCCESS
    except CancelledException:
        evaluation.status = models.Status.CANCELLED
        evaluation.append(stderr="Cancelled")
    except scoring.ScoringFailureException as exc:
        evaluation.status = models.Status.FAILURE
        evaluation.append(stderr=f"Error scoring: {exc}")
    except Exception as exc:  # pylint: disable=broad-except
        evaluation.status = models.Status.FAILURE
        evaluation.append(stderr=f"Execution failure: {exc}\n")
    finally:
        evaluation.save(update_fields=["status"])
        evaluation.cleanup_local_outputs(output_file_keys)

    return evaluation.status
