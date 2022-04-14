#! /bin/bash
#SBATCH --job-name=submtr
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --account=DMOBLEY_LAB
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=FAIL
#SBATCH --error=/data/homezvol0/osatom/sampl-app-extras/logs/start_remote_scheduler.sh.e
#SBATCH --output=/data/homezvol0/osatom/sampl-app-extras/logs/start_remote_scheduler.sh.o
#SBATCH --open-mode=append

echo "#################################################"
echo JOB ID: $SLURM_JOBID


echo "#################################################" 1>&2
echo JOB ID: $SLURM_JOBID 1>&2

. ~/.bashrc

echo CD into SAMPL-league/app
cd /data/homezvol0/osatom/SAMPL-league/app

echo Setting environment variables
source /data/homezvol0/osatom/sampl-app-extras/env/setlocal.sh

echo Running job_submitter.py
pipenv run python referee/job_submitter.py
