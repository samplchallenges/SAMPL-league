#!/bin/bash

if [[ ${ENV_TYPE} == 'production' ]]; then
    JOB_NAME='submtr'
elif [[ ${ENV_TYPE} == 'staging' ]]; then
    JOB_NAME='sbmtrstg'
else
    echo "Unknown ENV_TYPE: ${ENV_TYPE}" 1>&2
    exit 1
fi

if [[ ${SUBMITTER_QUEUE_PARTITION} == 'free' ]]; then
    sbatch \
        --job-name=$JOB_NAME \
        --partition=free \
        --output="${SAMPL_LOGS_ROOT}/staging-remote_scheduler.out" \
        --open-mode=append \
        "${SAMPL_ROOT}/app/remote_scheduler.slurm"
elif [[ ${SUBMITTER_QUEUE_PARTITION} == 'standard' ]]; then
    sbatch \
        --job-name=$JOB_NAME \
        --partition=standard \
        --account=DMOBLEY_LAB \
        --output="${SAMPL_LOGS_ROOT}/staging-remote_scheduler.out" \
        "${SAMPL_ROOT}/app/remote_scheduler.slurm"
else
    echo "Unknown SUBMITTER_QUEUE_PARTITION: ${SUBMITTER_QUEUE_PARTITION}" 1>&2
    exit 2
fi

