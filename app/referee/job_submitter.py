import logging
import subprocess
import time

import django
import yaml
from dask.distributed import Client
from dask_jobqueue import SLURMCluster

django.setup()
from django.conf import settings

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


def start_cluster(config_file, preload_file, worker_outfile, min_workers, max_workers):
    with open(config_file, encoding="utf-8") as f:
        config = yaml.safe_load(f)

    job_extra = [
        f"--output={worker_outfile}",
        "--open-mode=append",
    ]
    if settings.WORKER_QUEUE_PARTITION == "free":
        job_extra.append("--partition=free")
    elif settings.WORKER_QUEUE_PARTITION == "standard":
        job_extra.append("--partition=standard")
        job_extra.append("--account=DMOBLE_LAB")
    else:
        raise Exception(
            f"Unsupported WORKER_QUEUE_PARTITION {settings.WORKER_QUEUE_PARTITION}"
        )

    cluster = SLURMCluster(
        extra=[
            f"--preload {preload_file}",
        ],
        job_extra=job_extra,
        walltime=settings.WORKER_WALLTIME,
        **config["jobqueue"]["slurm"],
    )
    cluster.adapt(
        minimum_jobs=min_workers,
        maximum_jobs=max_workers,
    )
    return cluster


def set_submission_run_status(submission_run, status):
    submission_run.status = status
    submission_run.save(update_fields=["status"])


def check_for_submission_runs(start_time, client, check_interval, job_lifetime):
    n = 0
    logger.debug(
        "Checking for submissions every %d seconds over %d seconds",
        check_interval,
        job_lifetime,
    )
    while time.time() - start_time + (1.5 * check_interval) < job_lifetime:
        logger.debug("Checking for submission runs n=%d", n)
        for run in SubmissionRun.objects.filter(status=Status.PENDING_REMOTE):
            if run.is_public:
                logger.debug("Added run=%d", run.id)
                set_submission_run_status(run, Status.PENDING)
                rt.submit_submission_run(client, run)

            else:  # run is private
                for pair in SubmissionRunPair.objects.filter(private_run=run):
                    if pair.public_run.status == Status.SUCCESS:
                        set_submission_run_status(run, Status.PENDING)
                        rt.submit_submission_run(client, run)
                    elif pair.public_run.status in [
                        Status.CANCELLED,
                        Status.CANCEL_PENDING,
                        Status.FAILURE,
                    ]:
                        set_submission_run_status(run, Status.CANCELLED)
                    else:  # Status.PENDING_REMOTE or Status.RUNNING
                        pass

        time.sleep(check_interval)
        n += 1


def reset_unfinished_to_pending_submission():
    for status in [Status.PENDING, Status.RUNNING]:
        for submission_run in SubmissionRun.objects.filter(status=status):
            logger.debug(
                "Resetting PENDING/RUNNING submission_runs back to PENDING_REMOTE: %d",
                submission_run.id,
            )
            set_submission_run_status(submission_run, Status.PENDING_REMOTE)
            logger.debug("SubmissionRun status is now: %s", submission_run.status)
            for evaluation in submission_run.evaluation_set.all():
                if evaluation.status == Status.RUNNING:
                    logger.debug(
                        "   Resetting RUNNING back to PENDING: %d", evaluation.id
                    )
                    evaluation.status = Status.PENDING
                    evaluation.save(update_fields=["status"])
                    logger.debug("   Evaluation status is now: %s", evaluation.status)


def job_submitter_main():
    start_time = time.time()
    logger.info("Starting job_submitter.py at %s", time.ctime(start_time))
    logger.info("Container engine: %s", settings.CONTAINER_ENGINE)
    jobqueue_config_file = settings.SAMPL_ROOT / "app/referee/jobqueue.yaml"
    cluster = start_cluster(
        jobqueue_config_file,
        f"{settings.SAMPL_ROOT}/app/daskworkerinit.py",
        f"{settings.DASK_WORKER_LOGS_ROOT}/dask-worker-%j.out",
        settings.MINIMUM_WORKERS,
        settings.MAXIMUM_WORKERS,
    )
    logger.debug("Min: %d; Max: %d", settings.MINIMUM_WORKERS, settings.MAXIMUM_WORKERS)
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
