import subprocess
import smtplib
from datetime import datetime
from pytz import timezone

sender = 'from@fromdomain.com'
receivers = ['sampl-devops-aaaaexi2kkbdjwkxsqltscjpuy@mobleylab.slack.com']

date_format = '%Y-%m-%d %H:%M:%S %Z%z'
current_time = datetime.now(timezone('America/Los_Angeles'))
current_tstr = current_time.strftime(date_format)


message = f"""From: HPC3 Cron Job <hpc3cron@app.samplchallenges.org>
To: sampl-devops Slack Channel <sampl-devops-aaaaexi2kkbdjwkxsqltscjpuy@mobleylab.slack.com>
Subject: HPC3 ERROR: Scheduler is DOWN

The SLURMQueue Scheduler is DOWN
As of time: {current_tstr}
"""
def send_scheduler_down_email():
    try:
        smtpObj = smtplib.SMTP('localhost')
        smtpObj.sendmail(sender, receivers, message)         
        print("Successfully sent email")
    except SMTPException:
        print("Error: unable to send email")

def check_submitter_is_running():
    squeue_output = subprocess.check_output(["squeue", "-u", "osatom"]).decode('ascii')
    
    for line in squeue_output.split("\n"):
        if len(line) != 0:
            if line.split()[2] == "submtr":
                print("Job submitter is running or queued")
                return 0 
    
    send_scheduler_down_email()
    print("Job submitter is down")
    return 1

if __name__ == "__main__":
    check_submitter_is_running()
    
