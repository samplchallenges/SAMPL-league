#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=FAIL

echo "#################################################"
echo JOB ID: $SLURM_JOBID

. ~/.bashrc

echo CD into SAMPL-league/app
cd "${SAMPL_ROOT}/app"

echo Setting environment variables
source ${ENV_VAR_SCRIPT}

echo Running job_submitter.py
pipenv run python referee/job_submitter.py
