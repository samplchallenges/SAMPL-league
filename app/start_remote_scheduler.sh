#!/bin/bash

if [[ ${ENV_TYPE} == 'production' ]]; then
    JOB_NAME='submtr'
elif [[ ${ENV_TYPE} == 'staging' ]]; then
    JOB_NAME='sbmtrstg'
else
    echo "Unknown ENV_TYPE: ${ENV_TYPE}" 1>&2
    exit 1
fi

LOG_FILE="${SAMPL_LOGS_ROOT}/${ENV_TYPE}-remote_scheduler.out"

if [[ ${SUBMITTER_QUEUE_PARTITION} == 'free' ]]; then
    PARTITION="--partition=free"
elif [[ ${SUBMITTER_QUEUE_PARTITION} == 'standard' ]]; then
    PARTITION="--partition=standard --account=DMOBLEY_LAB"
else
    echo "Unknown SUBMITTER_QUEUE_PARTITION: ${SUBMITTER_QUEUE_PARTITION}" 1>&2
    exit 2
fi


sbatch \
    --job-name=$JOB_NAME \
    --time=$JOB_SUBMITTER_WALLTIME \
    $PARTITION \
    --output=$LOG_FILE \
    "${SAMPL_ROOT}/app/remote_scheduler.slurm"
