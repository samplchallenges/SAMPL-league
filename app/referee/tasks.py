import logging
import tempfile
import time
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
    public_run = submission.create_run(
        is_public=True, status=models.Status.PENDING_REMOTE
    )
    private_run = submission.create_run(
        is_public=False, status=models.Status.PENDING_REMOTE
    )
    submission.create_run_pair(public_run=public_run, private_run=private_run)


def run_and_score_submission(client, submission):
    """
    Runs public and private, plus scoring
    """
    delayed_conditional = dask.delayed(True)

    public_run = submission.create_run(is_public=True, status=models.Status.PENDING)
    private_run = submission.create_run(is_public=False, status=models.Status.PENDING)
    submission.create_run_pair(public_run=public_run, private_run=private_run)

    for submission_run in (public_run, private_run):
        delayed_conditional = _run(submission_run, delayed_conditional)

    if settings.VISUALIZE_DASK_GRAPH:
        delayed_conditional.visualize(filename="task_graph.svg")
    # return delayed_conditional.compute(scheduler="synchronous")

    future = client.submit(delayed_conditional.compute)  # pylint:disable=no-member
    logger.info("Future key: %s", future.key)

    dd.fire_and_forget(future)
    return future


def submit_submission_run(client, submission_run):
    delayed_conditional = dask.delayed(True)
    delayed_conditional = _run(submission_run, delayed_conditional)
    if settings.VISUALIZE_DASK_GRAPH:
        delayed_conditional.visualize(filename="task_graph.svg")
    future = client.submit(delayed_conditional.compute)  # pylint:disable=no-member
    logger.info("Future key: %s", future.key)

    dd.fire_and_forget(future)
    return future


@dask.delayed(pure=False)
def cache_containers(submission_run, delayed_conditional):

    if not delayed_conditional:
        return 1

    container = submission_run.submission.container
    aws_login_func = (
        utils.get_aws_credential_function(container.uri)
        if settings.LOGIN_TO_AWS
        else None
    )
    submission_run.append(stdout=f"{container.uri}\n")
    submission_run.append(stdout=f"{container.container_type}\n")

    if settings.CONTAINER_FILES_ROOT:
        container_save_path = (
            settings.CONTAINER_FILES_ROOT
            / f"{container.label.replace('/', '_')}_{container.tag}.sif"
            if container.tag
            else settings.CONTAINER_FILES_ROOT
            / f"{container.label.replace('/', '_')}.sif"
        )
    else:
        container_save_path = None

    pull_code, stdout, stderr = ever_given.wrapper.pull_container(
        container.uri,
        container.container_type,
        settings.CONTAINER_ENGINE,
        container_save_path,
        aws_login_func,
    )
    submission_run.append(stderr=stderr)
    submission_run.append(stdout=stdout)
    submission_run.append(stdout=f"PULL EXIT CODE: {pull_code} {type(pull_code)}\n")
    return pull_code


def _run(submission_run, delayed_conditional):

    pull_code = cache_containers(submission_run, delayed_conditional)

    if submission_run.submission.challenge.current_batch_group() is None:
        statuses = _run_evaluations(submission_run, pull_code, delayed_conditional)
    else:
        statuses = _run_batches(submission_run, pull_code, delayed_conditional)
    return check_and_score(submission_run.id, delayed_conditional, statuses)


def get_submission_run_status(statuses, submission_run_id):
    submission_run = models.SubmissionRun.objects.get(pk=submission_run_id)
    uniq_statuses = set(statuses)
    submission_run.append(stdout=f"Statuses {uniq_statuses}")
    if {models.Status.PENDING, models.Status.RUNNING} & uniq_statuses:
        submission_run.append(
            stderr=f"Evaluations should have all completed, but have statuses {statuses}!"
        )
        status = models.Status.FAILURE
    elif {models.Status.CANCELLED} == uniq_statuses or {
        models.Status.CANCELLED,
        models.Status.SUCCESS,
    } == uniq_statuses:
        status = models.Status.CANCELLED
    elif {models.Status.FAILURE, models.Status.CANCELLED} & uniq_statuses:
        status = models.Status.FAILURE
    elif {models.Status.SUCCESS} == uniq_statuses:
        status = models.Status.SUCCESS
    elif not uniq_statuses:
        submission_run.append(stderr="No sub statuses")
        status = models.Status.FAILURE
    else:
        submission_run.append(stderr=f"Sub statuses {uniq_statuses}")
        status = models.Status.FAILURE

    if (
        submission_run.status == models.Status.CANCEL_PENDING
        and status == models.Status.SUCCESS
    ):
        status = models.Status.CANCELLED

    return status


@dask.delayed(pure=False)  # pylint:disable=no-value-for-parameter
def check_and_score(submission_run_id, delayed_conditional, evaluation_statuses):
    submission_run = models.SubmissionRun.objects.get(pk=submission_run_id)
    if delayed_conditional:
        status = get_submission_run_status(evaluation_statuses, submission_run_id)
    else:
        status = models.Status.CANCELLED

    submission_run.status = status
    if status != models.Status.SUCCESS:
        submission_run.append(stderr=f"Submission run failed {status}")
        submission_run.save(update_fields=["status"])
        return False
    submission_run.append(stdout="Running check_and_score")
    scoring.score_submission_run(submission_run)
    return True


def _run_evaluations_or_batches(
    submission_run, pull_code, conditional, object_set, run_func
):
    statuses = []
    for obj in object_set.all():
        if obj.status == models.Status.PENDING:
            statuses.append(
                run_func(
                    submission_run.submission.id,
                    obj.id,
                    submission_run.id,
                    pull_code,
                    conditional=conditional,
                )
            )
        else:
            statuses.append(obj.status)
    return statuses


def _run_batches(submission_run, pull_code, conditional):
    return _run_evaluations_or_batches(
        submission_run,
        pull_code,
        conditional,
        submission_run.batchevaluation_set,
        run_batch_evaluation,
    )


def _run_evaluations(submission_run, pull_code, conditional):
    return _run_evaluations_or_batches(
        submission_run,
        pull_code,
        conditional,
        submission_run.evaluation_set,
        run_evaluation,
    )


def run_eval_or_batch(
    submission_id, cls, object_id, submission_run_id, pull_code, conditional
):
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

    if pull_code != 0:
        obj.status = models.Status.FAILURE
        obj.save(
            update_fields=[
                "status",
            ]
        )
        obj.append(stderr=f"Failed to pull container: {container.uri}")
        return models.Status.FAILURE

    output_file_keys = challenge.output_file_keys()
    all_output_keys = challenge.output_keys()

    if isinstance(
        obj, models.Evaluation
    ):  # or duck typing hasattr(obj, "input_element")
        kwargs, file_kwargs = values_helper.all_values(obj.input_element)
        is_batch = False
    else:
        kwargs, file_kwargs = values_helper.batch_values(obj.input_batch)
        is_batch = True

    obj.mark_started()
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

            if output_file_keys or is_batch:
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
                output_file_keys=all_output_keys if is_batch else output_file_keys,
                log_handler=cls.LogHandler(obj),
                cancel_requested_func=submission_run.check_cancel_requested,
                aws_login_func=aws_login_func,
            )

            for key, value in parsed_results:
                output_type = challenge.output_type(key)
                if not output_type:
                    obj.append(stderr=f"Ignoring key {key} with value {value}\n")
                    continue
                obj.clear_old_predictions(output_type)
                if is_batch:
                    models.Prediction.load_batch_output(
                        challenge, obj, output_type, value
                    )
                else:
                    models.Prediction.load_evaluation_output(
                        challenge, obj, output_type, value
                    )

        if is_batch:
            input_elements = list(obj.input_batch.elements())
        else:
            input_elements = [obj.input_element]

        obj.clear_old_scores()
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
def run_evaluation(
    submission_id, evaluation_id, submission_run_id, pull_code, conditional
):
    return run_eval_or_batch(
        submission_id,
        models.Evaluation,
        evaluation_id,
        submission_run_id,
        pull_code,
        conditional,
    )


@dask.delayed(pure=False)  # pylint:disable=no-value-for-parameter
def run_batch_evaluation(
    submission_id, batch_evaluation_id, submission_run_id, pull_code, conditional
):
    return run_eval_or_batch(
        submission_id,
        models.BatchEvaluation,
        batch_evaluation_id,
        submission_run_id,
        pull_code,
        conditional,
    )
