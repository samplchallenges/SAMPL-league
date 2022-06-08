import os
import time
from unittest.mock import Mock, patch

import dask.distributed as dd
import pytest
from django.contrib.contenttypes.models import ContentType
from django.core.management import call_command
from django.db import transaction

from core import models
from referee import job_submitter, tasks


def test_start_cluster():
    with patch("django.conf.settings.WORKER_QUEUE_PARTITION", "free"):
        config_file = f"{os.path.dirname(os.path.abspath(__file__))}/jobqueue_test.yaml"
        cluster = job_submitter.start_cluster(
            config_file, "SAMPL-league/app/daskworkerinit.py", "", 0, 2
        )
        job_script = cluster.job_script()
        assert "--partition=free" in job_script
        assert "--mem=4G" in job_script
        assert "--preload SAMPL-league/app/daskworkerinit.py" in job_script
        assert "--cpus-per-task=1" in job_script


@pytest.mark.django_db(transaction=True)
def test_reset_unfinished_to_pending_submission():
    processes = True
    if processes:
        transaction.commit()
        call_command("migrate", "core", "zero", interactive=False)
        call_command("migrate", "core", interactive=False)
        call_command("sample_data")
        transaction.commit()
    else:
        call_command("sample_data")

    submission = models.Submission.objects.first()
    tasks.enqueue_submission(submission)

    # set up SubmissionRun and Evaluation statuses to be modified by
    # job_submitter.reset_unfinished_to_pending_submission()
    submission_run1 = models.SubmissionRun.objects.get(pk=2)
    submission_run1.status = models.Status.RUNNING
    submission_run1.save(update_fields=["status"])
    for evaluation in submission_run1.evaluation_set.all():
        evaluation.status = models.Status.RUNNING
        evaluation.save(update_fields=["status"])

    submission_run2 = models.SubmissionRun.objects.get(pk=3)
    submission_run2.status = models.Status.PENDING
    submission_run2.save(update_fields=["status"])

    job_submitter.reset_unfinished_to_pending_submission()

    # check that SubmissionRun and Evaluation statuses were modified
    submission_run1 = models.SubmissionRun.objects.get(pk=2)
    submission_run2 = models.SubmissionRun.objects.get(pk=3)
    assert submission_run1.status == models.Status.PENDING_REMOTE
    for evaluation in submission_run1.evaluation_set.all():
        assert evaluation.status == models.Status.PENDING

    assert submission_run2.status == models.Status.PENDING_REMOTE


@pytest.mark.parametrize(
    ["container_engine"],
    [["docker"], ["singularity"]],
)
@pytest.mark.django_db(transaction=True)
def test_check_for_submission_runs(client, container_engine):
    with patch("django.conf.settings.CONTAINER_ENGINE", container_engine):
        processes = True
        if processes:
            transaction.commit()
            call_command("migrate", "core", "zero", interactive=False)
            call_command("migrate", "core", interactive=False)
            call_command("sample_data")
            transaction.commit()
        else:
            call_command("sample_data")

        submission = models.Submission.objects.first()
        tasks.enqueue_submission(submission)

        os.environ["CONTAINER_ENGINE"] = container_engine
        preload_file = "daskworkerinit_tst.py"
        cluster = dd.LocalCluster(n_workers=0, preload=(preload_file,))
        dask_client = dd.Client(cluster)

        job_submitter.check_for_submission_runs(time.time(), dask_client, 1, 3)

        # time.sleep(120)
        # public run
        submission_run = models.SubmissionRun.objects.get(pk=2)
        assert submission_run.status == models.Status.PENDING
        # private run
        submission_run = models.SubmissionRun.objects.get(pk=3)
        assert submission_run.status == models.Status.PENDING_REMOTE
        # submission_run = models.SubmissionRun.objects.first()
        # assert submission_run.status == models.Status.SUCCESS


def test_resubmit_check_for_submission_runs_job():
    with patch("subprocess.check_output", mock_check_output):
        job_info = job_submitter.resubmit_check_for_submission_runs_job(
            "job_submission_script.sh"
        )
        assert job_info == "Submitted batch job 11163505\n"


def mock_check_output(*args, **kwargs):
    return b"Submitted batch job 11163505\n"
