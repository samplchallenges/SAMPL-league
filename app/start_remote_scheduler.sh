#! /bin/bash
#SBATCH --job-name=test-pipenv
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --account=osatom
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-type=FAIL
#SBATCH --error=/data/homezvol0/osatom/sampl-app-extras/logs/start_remote_scheduler.sh.e
#SBATCH --output=/data/homezvol0/osatom/sampl-app-extras/logs/start_remote_scheduler.sh.o
#SBATCH --open-mode=append

echo JOB ID: $SLURM_JOBID
echo JOB ID: $SLURM_JOBID 1>&2

. ~/.bashrc

cd /data/homezvol0/osatom/SAMPL-league/app
echo cd SAMPL-league/app

echo activating environment
pipenv shell

echo setting environment variables
source /data/homezvol0/osatom/sampl-app-extras/env/setlocal.sh

pip freeze > envout

echo running job_submitter.py
pipenv run python referee/job_submitter.py >> outfile 2>&1

exit 0
