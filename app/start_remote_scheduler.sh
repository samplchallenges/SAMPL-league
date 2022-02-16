#! /bin/bash
#SBATCH --job-name=test-pipenv
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --account=DMOBLEY_LAB
#SBATCH --error=pop.e
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=FAIL

. ~/.bashrc

cd $SLURM_SUBMIT_DIR
pipenv shell

source setlocal.sh

pip freeze > envout

#django-admin runserver >> django-server.log 2>&1 &

pipenv run dask-scheduler --host 0.0.0.0:8786  >> /data/homezvol0/osatom/sampl-app-extras/logs/dask-scheduler.log 2>&1 &

pipenv run dask-worker --nthreads 1 --nprocs 1 --preload daskworkerinit.py 127.0.0.1:8786 >> /data/homezvol0/osatom/sampl-app-extras/logs/dask-worker.log 2>&1 &

pipenv run python referee/job_submitter.py >> outfile 2>&1

exit 0
