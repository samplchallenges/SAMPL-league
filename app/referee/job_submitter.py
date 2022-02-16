import os
import sys
import subprocess
import time

from dask.distributed import Client

import django
django.setup()
from django.conf import settings

from core.models import Status, SubmissionRun
import referee.tasks as rt


JOB_SUBMITTER_LIFETIME = 3600 # seconds
CHECK_INTERVAL = 120 # seconds



def _resubmit_check_for_submission_runs_job():

    os.system("qsub /data/homezvol0/osatom/SAMPL-league/app/start_remote_scheduler.sh")

    print("resubmitted job to queue")
    sys.stdout.flush()

def inc(x):
    return x + 1

def check_for_submission_runs(): 
    time_started = time.time()
    num_intervals = 1

    client = Client(settings.DASK_SCHEDULER_URL)
    #x_val = 10
    #print("Client created")
    #print("\tinit: x_val = {} should be 10".format(x_val))
    #sys.stdout.flush()

    # Any jobs that status.PENDING or status.RUNNING at this point are no longer queued, and need
    # to be re-queued at start up

    for run in SubmissionRun.objects.filter(status=Status.RUNNING):
        # submit jobs
        run.status = Status.PENDING
        run.save(update_fields=["status"])

    for run in SubmissionRun.objects.filter(status=Status.PENDING):
        # submit jobs
        pass

    while time.time() - time_started + CHECK_INTERVAL < JOB_SUBMITTER_LIFETIME:
        
        #print("check interval {}".format(num_intervals))
        #print("\ttime:", time.time() - time_started)
        #x = client.submit(inc, x_val)
        #x_val = x.result()
        #print("\tresult: x_val = {} should be 10 + {}".format(x_val, num_intervals))
        #sys.stdout.flush()

        #import urllib.request
        #external_ip = urllib.request.urlopen('https://ident.me').read().decode('utf8')
        #print("\texternal ip:", external_ip)


        print("\tprinting pending submissions")
        for run in SubmissionRun.objects.filter(status=Status.CANCELLED):
            if run.id != 523:
                continue
            print("\t\t", run)
            run.status = Status.PENDING
            run.save(update_fields=["status"])
            rt.submit_submission_run(client, run)
        print("\tfinished printing pending submissions")
        sys.stdout.flush()
        time.sleep(CHECK_INTERVAL)
        sys.stdout.flush()
        num_intervals += 1
    

    _resubmit_check_for_submission_runs_job() 
    print("exiting this execution")


if __name__ == "__main__":
    check_for_submission_runs()
