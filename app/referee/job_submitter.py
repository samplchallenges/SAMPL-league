import logging
import subprocess
import time

import django
import yaml
from dask.distributed import Client
from dask_jobqueue import SLURMCluster

django.setup()
from django.conf import settings
from django.db import transaction

import referee.tasks as rt
from core.models import Status, SubmissionRun, SubmissionRunPair

logger = logging.getLogger("job_submitter")
logger.setLevel(logging.DEBUG)


def resubmit_check_for_submission_runs_job(scheduler_submission_script):
    jobid_sub_bytes = subprocess.check_output(
        ["sh", scheduler_submission_script],
    )
    # jobid_sub_bytes is a bytes string that looks like b'Submitted batch job 11163505\n'
    return jobid_sub_bytes.decode("utf-8")


def start_cluster(config_file, preload_file, worker_outfile):
    with open(config_file, encoding="utf-8") as f:
        config = yaml.safe_load(f)

    job_extra = [
        f"--output={worker_outfile}",
        "--open-mode=append",
    ]

    cluster = SLURMCluster(
        extra=[
            f"--preload {preload_file}",
        ],
        job_extra=job_extra,
        **config["jobqueue"]["slurm"]["cluster_settings"],
    )
    logger.debug(
        "Min: %d; Max: %d",
        config["jobqueue"]["slurm"]["adapt_settings"]["minimum_jobs"],
        config["jobqueue"]["slurm"]["adapt_settings"]["maximum_jobs"],
    )
    cluster.adapt(**config["jobqueue"]["slurm"]["adapt_settings"])
    return cluster


def set_submission_run_status(submission_run, status):
    submission_run.status = status
    submission_run.save(update_fields=["status"])


def job_submitter_alive(start_time, check_interval, job_lifetime):
    return time.time() - start_time + (3 * check_interval) < job_lifetime


def check_for_submission_runs(start_time, client, check_interval, job_lifetime):
    n = 0
    logger.debug(
        "Checking for submissions every %d seconds over %d seconds",
        check_interval,
        job_lifetime,
    )

    while job_submitter_alive(start_time, check_interval, job_lifetime):
        logger.debug("Checking for submission runs n=%d", n)
        with transaction.atomic():
            submission_run_query = (
                SubmissionRun.objects.select_for_update()
                .filter(status=Status.PENDING_REMOTE)
                .all()
            )
            for run in submission_run_query:
                if run.is_public:
                    logger.debug("Added run=%d", run.id)
                    run.update_status(Status.PENDING)
                    rt.submit_submission_run(client, run)

                else:  # run is private
                    for pair in SubmissionRunPair.objects.filter(private_run=run):
                        if pair.public_run.status == Status.SUCCESS:
                            run.update_status(run, Status.PENDING)
                            rt.submit_submission_run(client, run)
                        elif pair.public_run.status in [
                            Status.CANCELLED,
                            Status.CANCEL_PENDING,
                            Status.FAILURE,
                        ]:
                            run.update_status(Status.CANCELLED)
                        else:  # Status.PENDING_REMOTE or Status.RUNNING
                            pass
        if job_submitter_alive(start_time, check_interval, job_lifetime):
            time.sleep(check_interval)
        n += 1


def reset_unfinished_to_pending_submission():
    for status in [Status.PENDING, Status.RUNNING]:
        with transaction.atomic():
            submission_run_query = (
                SubmissionRun.objects.select_for_update().filter(status=status).all()
            )
            for submission_run in submission_run_query:
                logger.debug(
                    "Resetting PENDING/RUNNING submission_runs back to PENDING_REMOTE: %d",
                    submission_run.id,
                )
                submission_run.update_status(Status.PENDING_REMOTE)
                logger.debug("SubmissionRun status is now: %s", submission_run.status)

                evaluation_query = (
                    submission_run.evaluation_set.select_for_update().all()
                )
                for evaluation in evaluation_query:
                    if evaluation.status == Status.RUNNING:
                        logger.debug(
                            "   Resetting RUNNING back to PENDING: %d", evaluation.id
                        )
                        evaluation.update_status(Status.PENDING)
                        logger.debug(
                            "   Evaluation status is now: %s", evaluation.status
                        )


def job_submitter_main():
    start_time = time.time()
    logger.info("Starting job_submitter.py at %s", time.ctime(start_time))
    logger.info("Job Lifetime: %s", settings.JOB_SUBMITTER_LIFETIME)
    logger.info("Container engine: %s", settings.CONTAINER_ENGINE)
    jobqueue_config_file = settings.JOBQUEUE_CONFIG_FILE
    cluster = start_cluster(
        jobqueue_config_file,
        f"{settings.SAMPL_ROOT}/app/daskworkerinit.py",
        f"{settings.DASK_WORKER_LOGS_ROOT}/dask-worker-%j.out",
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

    scheduler_submission_script = settings.SAMPL_ROOT / "app/start_remote_scheduler.sh"
    jobinfo = resubmit_check_for_submission_runs_job(scheduler_submission_script)
    logger.info("Resubmitting start_remote_scheduler.sh - JOB INFO: %s", jobinfo)


if __name__ == "__main__":
    job_submitter_main()
