#!/bin/bash

sbatch \
    --output="${SAMPL_LOGS_ROOT}/remote_scheduler.out" \
    --open-mode="append" \
    "${SAMPL_ROOT}/app/remote_scheduler.slurm"
