import logging
import os
import subprocess
import sys
import time

import django
import yaml
from dask.distributed import Client
from dask_jobqueue import SLURMCluster

django.setup()
from django.conf import settings

import referee.tasks as rt
from core.models import Status, SubmissionRun

logger = logging.getLogger("job_submitter")
logger.setLevel(logging.DEBUG)


def resubmit_check_for_submission_runs_job():
    scheduler_submission_script = settings.SAMPL_ROOT / "app/start_remote_scheduler.sh"
    jobid_sub_bytes = subprocess.check_output(
        f"sbatch {scheduler_submission_script}", shell=True
    )
    # jobid_sub_bytes is a bytes string that looks like b'Submitted batch job 11163505\n'
    return str(jobid_sub_bytes)


def start_cluster(jobqueue_config_file):
    with open(jobqueue_config_file) as f:
        config = yaml.safe_load(f)
    cluster = SLURMCluster(**config["jobqueue"]["slurm"])
    return cluster


def check_for_submission_runs(start_time, client, check_interval, job_lifetime):
    n = 0
    logger.debug("Checking for submissions every %d seconds over %d seconds")
    while time.time() - start_time + (1.5 * check_interval) < job_lifetime:
        logger.debug("Checking for submission runs n=%d", n)
        for run in SubmissionRun.objects.filter(status=Status.CANCELLED):
            if run.id > 760:
                logger.debug("Added run=%d", run.id)
                run.status = Status.PENDING
                run.save(update_fields=["status"])

                for evaluation in run.evaluation_set.all():
                    evaluation.status = Status.PENDING
                    evaluation.save(update_fields=["status"])
                rt.submit_submission_run(client, run)
        time.sleep(check_interval)
        n += 1


def reset_unfinished_to_pending_submission():
    for status in [Status.PENDING, Status.RUNNING]:
        for submission_run in SubmissionRun.objects.filter(status=status):
            logger.debug(
                "Resetting PENDING/RUNNING submission_runs back to PENDING_REMOTE: %d",
                submission_run.id,
            )
            submission_run.status = Status.PENDING_REMOTE
            submission_run.save(update_fields=["status"])
            for evaluation in submission_run.evaluation_set.all():
                if evaluation.status == Status.RUNNING:
                    logger.debug(
                        "   Resetting RUNNING back to PENDING: %d", evaluation.id
                    )
                    evaluation.status = Status.PENDING
                    evaluation.save(update_fields=["status"])


if __name__ == "__main__":
    start_time = time.time()
    logger.info("Starting job_submitter.py at %s", time.ctime(start_time))
    logger.info("Container engine: %s", settings.CONTAINER_ENGINE)

    jobqueue_config_file = settings.SAMPL_ROOT / "app/referee/jobqueue.yaml"

    cluster = start_cluster(jobqueue_config_file)

    logger.debug("Min: %d; Max: %d", settings.MINIMUM_WORKERS, settings.MAXIMUM_WORKERS)
    cluster.adapt(
        minimum_jobs=settings.MINIMUM_WORKERS, maximum_jobs=settings.MAXIMUM_WORKERS
    )

    client = Client(cluster)

    check_for_submission_runs(
        start_time,
        client,
        settings.CHECK_INTERVAL,
        settings.JOB_SUBMITTER_LIFETIME,
    )

    # shutdown scheduler and connected workers
    client.shutdown()
    logger.info("Shut down Dask scheduler and workers")

    # switch the status of any jobs that are currently running back to
    # pending submission to queue so next instance of scheduler will requeue them
    reset_unfinished_to_pending_submission()

    jobinfo = resubmit_check_for_submission_runs_job()
    logger.info("Resubmitting start_remote_scheduler.sh - JOB INFO: %s", jobinfo)
