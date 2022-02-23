import os
import subprocess
import sys
import time
import logging
import yaml

from dask_jobqueue import SLURMCluster
from dask.distributed import Client

import django
django.setup()
from django.conf import settings

import referee.tasks as rt
from core.models import Status, SubmissionRun

logger = logging.getLogger("job_submitter")
logger.setLevel(logging.DEBUG)


def _resubmit_check_for_submission_runs_job():
    scheduler_submission_script = settings.SAMPL_ROOT / "app/start_remote_scheduler.sh"
    jobid = subprocess.check_output(
                f"sbatch {scheduler_submission_script}",
                shell = True
            )
    return jobid



def start_cluster(jobqueue_config_file):
    with open(jobqueue_config_file, 'r') as f:
        config = yaml.safe_load(f)
    cluster = SLURMCluster(**config['jobqueue']['slurm'])
    return cluster


def check_for_submission_runs(start_time, client, check_interval, job_lifetime):
    n = 0
    logger.debug(f"Checking for submissions every {check_interval} seconds over {job_lifetime} seconds")
    while time.time() - start_time + (1.5*check_interval) < job_lifetime:
        logger.debug(f"Checking for submission runs n={n}")
        for run in SubmissionRun.objects.filter(status=Status.CANCELLED):
            if run.id < 549 or run.id % 2 == 0:
                continue
            logger.debug(f"Added run={run.id}")
            run.status = Status.PENDING
            run.save(update_fields=["status"])
            for evaluation in run.evaluation_set.all():
                if evaluation.status == Status.CANCELLED:
                    evaluation.status = Status.PENDING
                    evaluation.save(update_fields=["status"])
            rt.submit_submission_run(client, run)
        time.sleep(check_interval)
        n += 1


def reset_unfinished_to_pending_submission():
    for submission_run in SubmissionRun.objects.filter(status=Status.PENDING):
        logger.debug(f"Resetting PENDING back to PENDING_REMOTE: {submission_run.id}")
        submission_run.status = Status.CANCELLED #Status.PENDING_REMOTE
        submission_run.save(update_fields=["status"])
        for evaluation in submission_run.evaluation_set.all():
            if evaluation.status == Status.RUNNING:
                logger.debug(f"   Resetting RUNNING back to PENDING: {evaluation.id}") 
                evaluation.status = Status.PENDING
                evaluation.save(update_fields=["status"])

if __name__ == "__main__":
    start_time = time.time()
    logger.info(f"Starting job_submitter.py at {time.ctime(start_time)}")

    jobqueue_config_file = settings.SAMPL_ROOT / "app/referee/jobqueue.yaml"
    
    cluster = start_cluster(jobqueue_config_file)

    logger.debug(f"Min: {settings.MINIMUM_WORKERS}; Max: {settings.MAXIMUM_WORKERS}") 
    cluster.adapt(
        minimum_jobs=settings.MINIMUM_WORKERS, 
        maximum_jobs=settings.MAXIMUM_WORKERS
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

    jobid = _resubmit_check_for_submission_runs_job()
    logger.info(f"Resubmitting start_remote_scheduler.sh - JOB ID: {jobid}")

