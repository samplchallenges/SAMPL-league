import time

import click

from core.models import Status, SubmissionRun
from referee.tasks import submit_submission_run

JOB_SUBMITTER_LIFETIME = 600


@click.command()
def check_for_submission_runs():
    time_started = time.now()
    while time.now() - time_started < JOB_SUBMITTER_LIFETIME:
        for run in SubmissionRun.objects.filter(status=Status.PENDING_REMOTE):
            run.status = Status.PENDING
            submit_submission_run(run)
            run.save(update_fields=["status"])
        time.sleep(CHECK_INTERVAL)

    resubmit_myself()
