import logging
import tempfile
from pathlib import Path

import dask
import dask.distributed as dd
import ever_given.wrapper
from django.conf import settings

from core import filecache, models

from . import scoring

logger = logging.getLogger(__name__)


def run_and_score_submission(client, submission):
    """
    Runs public and private, plus scoring
    """
    challenge = submission.challenge
    delayed_conditional = dask.delayed(True)
    for is_public in (True, False):
        element_ids = challenge.inputelement_set.filter(
            is_public=is_public
        ).values_list("id", flat=True)
        run_id, evaluation_statuses = build_submission_run(
            submission.pk, element_ids, delayed_conditional, is_public=is_public
        )
        delayed_conditional = check_and_score(
            run_id, delayed_conditional, evaluation_statuses
        )

    if settings.VISUALIZE_DASK_GRAPH:
        delayed_conditional.visualize(filename="task_graph.svg")

    future = client.submit(delayed_conditional.compute)  # pylint:disable=no-member
    logger.info("Future key: %s", future.key)

    dd.fire_and_forget(future)
    return future


@dask.delayed(pure=False)  # pylint:disable=no-value-for-parameter
def check_and_score(submission_run_id, delayed_conditional, evaluation_statuses):
    uniq_statuses = set(evaluation_statuses)
    if not delayed_conditional:
        status = models.Status.CANCELLED
    elif {models.Status.PENDING, models.Status.RUNNING} & uniq_statuses:
        logger.error(
            "Evaluations should have all completed, but have statuses %s!",
            evaluation_statuses,
        )
        status = models.Status.FAILURE
    elif {models.Status.CANCELLED} == uniq_statuses:
        status = models.Status.CANCELLED
    elif {models.Status.FAILURE, models.Status.CANCELLED} & uniq_statuses:
        status = models.Status.FAILURE
    else:
        status = models.Status.SUCCESS

    submission_run = models.SubmissionRun.objects.get(pk=submission_run_id)
    submission_run.status = status
    submission_run.save()
    if status != models.Status.SUCCESS:
        logger.warning(
            "Submission run failed (public? %s), %s: %s",
            submission_run.is_public,
            submission_run,
            status,
        )
        return False
    logger.info(
        "Running check_and_score %s public? %s",
        submission_run_id,
        submission_run.is_public,
    )
    scoring.score_submission(submission_run.submission.pk, submission_run_id)
    return True


def create_submission_run(submission_id, *, is_public):
    submission = models.Submission.objects.get(pk=submission_id)
    container = submission.container
    if not container.digest:
        container.digest = "nodigest"
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


def build_submission_run(submission_id, element_ids, conditional, is_public=True):

    submission_run_id = create_submission_run(submission_id, is_public=is_public)

    evaluations = [
        models.Evaluation.objects.create(
            input_element_id=element_id, submission_run_id=submission_run_id
        )
        for element_id in element_ids
    ]

    evaluation_statuses = [
        run_evaluation(
            submission_id,
            evaluation.id,
            submission_run_id,
            conditional=conditional,
        )
        for evaluation in evaluations
    ]
    return (submission_run_id, evaluation_statuses)


@dask.delayed(pure=False)  # pylint:disable=no-value-for-parameter
def run_evaluation(submission_id, evaluation_id, submission_run_id, conditional):
    if not conditional:
        models.Evaluation.objects.filter(pk=evaluation_id).update(
            status=models.Status.CANCELLED
        )
        return models.Status.CANCELLED
    submission = models.Submission.objects.get(pk=submission_id)
    container = submission.container
    challenge = submission.challenge
    submission_run = submission.submissionrun_set.get(pk=submission_run_id)
    evaluation_score_types = {
        score_type.key: score_type
        for score_type in challenge.scoretype_set.filter(
            level=models.ScoreType.Level.EVALUATION
        )
    }

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
            )

            for key, value in parsed_results:
                output_type = challenge.output_type(key)
                if output_type:
                    prediction = models.Prediction.load_output(
                        challenge, evaluation, output_type, value
                    )
                    evaluation.append(stdout=f"{prediction.__dict__}\n")
                    prediction.save()
                else:
                    evaluation.append(stderr=f"Ignoring key {key} with value {value}\n")

        try:
            scoring.score_evaluation(
                challenge.scoremaker.container,
                evaluation,
                evaluation_score_types,
            )
        except Exception as exc:  # pylint: disable=broad-except
            evaluation.append(stderr=f"Error scoring\n{exc}")
            evaluation.save(update_fields=["log_stderr"])
            raise
        evaluation.status = models.Status.SUCCESS
        return evaluation.status
    except Exception as exc:  # pylint: disable=broad-except
        evaluation.status = models.Status.FAILURE
        evaluation.append(stderr=f"Execution failure: {exc}\n")
        evaluation.save(update_fields=["log_stderr"])
        return evaluation.status
    finally:
        evaluation.save()
        for prediction in evaluation.prediction_set.filter(
            value_type__key__in=output_file_keys
        ):
            filecache.delete_local_cache(prediction.value)
