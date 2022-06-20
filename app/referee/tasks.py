import logging
import tempfile
from pathlib import Path

import dask
import dask.distributed as dd
import ever_given.wrapper
from django.conf import settings
from ever_given.log_processing import CancelledException

from core import models, values_helper

from . import scoring, utils

logger = logging.getLogger(__name__)


def enqueue_submission(submission):
    """
    Runs on the webapp in an environment that can't talk to the scheduler.
    Records in the database that we want the job submitter to call
    submit_submission_run
    """
    for is_public in (True, False):
        submission.create_run(is_public=is_public, status=models.Status.PENDING_REMOTE)


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


def submit_submission_run(client, submission_run):
    delayed_conditional = dask.delayed(True)
    delayed_conditional = _run(submission_run, delayed_conditional)
    future = client.submit(delayed_conditional.compute)  # pylint:disable=no-member
    logger.info("Future key: %s", future.key)

    dd.fire_and_forget(future)
    return future


def _trigger_submission_run(submission, delayed_conditional, *, is_public):
    submission_run = submission.create_run(
        is_public=is_public, status=models.Status.PENDING
    )
    return _run(submission_run, delayed_conditional)


def _run(submission_run, delayed_conditional):
    if submission_run.submission.challenge.current_batch_group() is None:
        statuses = _run_evaluations(submission_run, delayed_conditional)
    else:
        statuses = _run_batches(submission_run, delayed_conditional)
    return check_and_score(submission_run.id, delayed_conditional, statuses)


@dask.delayed(pure=False)  # pylint:disable=no-value-for-parameter
def check_and_score(submission_run_id, delayed_conditional, statuses):
    submission_run = models.SubmissionRun.objects.get(pk=submission_run_id)
    uniq_statuses = set(statuses)
    if not delayed_conditional:
        status = models.Status.CANCELLED
    elif {models.Status.PENDING, models.Status.RUNNING} & uniq_statuses:
        submission_run.append(
            stderr=f"Evaluations should have all completed, but have statuses {statuses}!"
        )
        status = models.Status.FAILURE
    elif {models.Status.CANCELLED} == uniq_statuses:
        status = models.Status.CANCELLED
    elif {models.Status.FAILURE, models.Status.CANCELLED} & uniq_statuses:
        status = models.Status.FAILURE
    elif {models.Status.SUCCESS} == uniq_statuses:
        status = models.Status.SUCCESS
    else:
        submission_run.append(f"Unrecognized statuses: {statuses}")
        status = models.Status.FAILURE
    submission_run.status = status
    if status != models.Status.SUCCESS:
        submission_run.append(stderr=f"Submission run failed {status}")
        submission_run.save(update_fields=["status"])
        return False
    submission_run.append(stdout="Running check_and_score")
    scoring.score_submission_run(submission_run)
    return True


def _run_evaluations_or_batches(submission_run, conditional, object_set, run_func):
    statuses = []
    for obj in object_set.all():
        if obj.status == models.Status.PENDING:
            statuses.append(
                run_func(
                    submission_run.submission.id,
                    obj.id,
                    submission_run.id,
                    conditional=conditional,
                )
            )
        else:
            statuses.append(obj.status)
    return statuses


def _run_batches(submission_run, conditional):
    return _run_evaluations_or_batches(
        submission_run,
        conditional,
        submission_run.batchevaluation_set,
        run_batch_evaluation,
    )


def _run_evaluations(submission_run, conditional):
    return _run_evaluations_or_batches(
        submission_run, conditional, submission_run.evaluation_set, run_evaluation
    )


def run_eval_or_batch(submission_id, cls, object_id, submission_run_id, conditional):
    submission = models.Submission.objects.get(pk=submission_id)
    container = submission.container
    challenge = submission.challenge
    submission_run = submission.submissionrun_set.get(pk=submission_run_id)
    if not conditional or submission_run.check_cancel_requested():
        cls.objects.filter(pk=object_id).update(status=models.Status.CANCELLED)
        return models.Status.CANCELLED
    if submission_run.status == models.Status.PENDING:
        submission_run.status = models.Status.RUNNING
        submission_run.save(update_fields=["status"])

    evaluation_score_types = challenge.score_types[models.ScoreType.Level.EVALUATION]
    obj = cls.objects.filter(submission_run_id=submission_run_id).get(pk=object_id)

    if obj.status not in {models.Status.PENDING, models.Status.RUNNING}:
        return obj.status

    output_file_keys = challenge.output_file_keys()

    if isinstance(
        obj, models.Evaluation
    ):  # or duck typing hasattr(obj, "input_element")
        kwargs, file_kwargs = values_helper.all_values(obj.input_element)
        is_batch = False
    else:
        kwargs, file_kwargs = values_helper.batch_values(obj)
        is_batch = True

    obj.mark_started(kwargs, file_kwargs)
    kwargs.update(container.custom_args())
    file_kwargs.update(container.custom_file_args())
    aws_login_func = (
        utils.get_aws_credential_function(container.uri)
        if settings.LOGIN_TO_AWS
        else None
    )
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            dirpath = Path(str(tmpdir))
            output_dir = None

            if output_file_keys:
                output_dir = dirpath / "output"
                output_dir.mkdir()

            parsed_results = ever_given.wrapper.run(
                container.uri,
                "--batch" if is_batch else "",
                container_type=container.container_type,
                engine_name=settings.CONTAINER_ENGINE,
                kwargs=kwargs,
                file_kwargs=file_kwargs,
                output_dir=output_dir,
                output_file_keys=output_file_keys,
                log_handler=cls.LogHandler(obj),
                cancel_requested_func=submission_run.check_cancel_requested,
                aws_login_func=aws_login_func,
            )

            for key, value in parsed_results:
                output_type = challenge.output_type(key)
                if output_type:
                    if is_batch:
                        models.Prediction.load_batch_output(
                            challenge, obj, output_type, value
                        )
                    else:
                        models.Prediction.load_evaluation_output(
                            challenge, obj, output_type, value
                        )
                else:
                    obj.append(stderr=f"Ignoring key {key} with value {value}\n")

        if is_batch:
            input_elements = []
        else:
            input_elements = [obj.input_element]
        for input_element in input_elements:
            # TODO: batch up scoring?
            for log_message in scoring.score_element(
                challenge.scoremaker.container,
                obj,
                submission_run,
                input_element,
                evaluation_score_types,
            ):
                obj.append(log_message)

        obj.status = models.Status.SUCCESS
    except CancelledException:
        obj.status = models.Status.CANCELLED
        obj.append(stderr="Cancelled")
    except scoring.ScoringFailureException as exc:
        obj.status = models.Status.FAILURE
        obj.append(stderr=f"Error scoring: {exc}")
    except Exception as exc:  # pylint: disable=broad-except
        obj.status = models.Status.FAILURE
        obj.append(stderr=f"Execution failure: {exc}\n")
    finally:
        obj.save(update_fields=["status"])
        obj.cleanup_local_outputs(output_file_keys)

    return obj.status


@dask.delayed(pure=False)  # pylint:disable=no-value-for-parameter
def run_evaluation(submission_id, evaluation_id, submission_run_id, conditional):
    return run_eval_or_batch(
        submission_id,
        models.Evaluation,
        evaluation_id,
        submission_run_id,
        conditional,
    )


@dask.delayed(pure=False)  # pylint:disable=no-value-for-parameter
def run_batch_evaluation(
    submission_id, batch_evaluation_id, submission_run_id, conditional
):
    return run_eval_or_batch(
        submission_id,
        models.BatchEvaluation,
        batch_evaluation_id,
        submission_run_id,
        conditional,
    )
