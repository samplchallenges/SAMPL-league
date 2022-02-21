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

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def _resubmit_check_for_submission_runs_job():
    scheduler_submission_script = settings.BASE_DIR / "SAMPL-league/app/start_remote_scheduler.sh"
    jobid = subprocess.check_output(f"sbatch {scheduler_submission_script}")
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
            
            rt.submit_submission_run(client, run)
        time.sleep(check_interval)
        n += 1


def reset_unfinished_to_pending_submission():
    for run in SubmissionRun.objects.filter(status=Status.PENDING):
        if run.id != 545 and run.id != 547:
            continue
        logger.debug(f"Resetting back to CANCELLED: {run.id}")
        run.status = Status.CANCELLED
        run.save(update_fields=["status"])
        rt.submit_submission_run(client, run)
    

if __name__ == "__main__":
    start_time = time.time()
    logger.info(f"Starting job_submitter.py at {start_time}")

    jobqueue_config_file = settings.BASE_DIR / "SAMPL-league/app/referee/jobqueue.yaml"
    cluster = start_cluster(jobqueue_config_file)
    
    logger.debug(f"Min: {settings.MINIMUM_WORKERS}; Max: {settings.MAXIMUM_WORKERS}") 
    cluster.adapt(
        minimum=settings.MINIMUM_WORKERS, 
        maximum=settings.MAXIMUM_WORKERS
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

