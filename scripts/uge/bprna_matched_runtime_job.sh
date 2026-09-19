#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=04:00:00
#$ -t 1-27
#$ -N bprna_matched_runtime
#$ -o logs/bprna_matched_runtime_$JOB_ID.$TASK_ID.out
#$ -e logs/bprna_matched_runtime_$JOB_ID.$TASK_ID.err

set -euo pipefail

REPO_ROOT="/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper"
TASK_ID="${SGE_TASK_ID:?SGE_TASK_ID is required}"
OUTDIR="${REPO_ROOT}/benchmark/results/bprna_matched_runtime/${JOB_ID:-manual}"
CHUNK_FILE="${REPO_ROOT}/benchmark/data/bprna_chunks/chunk_$(printf '%02d' "$TASK_ID").fa"
LINRACC_BIN="${REPO_ROOT}/vendor/LinearRaccess/build-raccess/LinRacc"
RACCESS_BIN="${REPO_ROOT}/vendor/raccess/src/raccess/run_raccess"

mkdir -p "$OUTDIR" "${REPO_ROOT}/benchmark/logs"

cd "$REPO_ROOT"
exec python3 "${REPO_ROOT}/benchmark/scripts/run_bprna_chunk_matched.py" \
    --chunk "$CHUNK_FILE" \
    --outdir "$OUTDIR" \
    --linracc "$LINRACC_BIN" \
    --raccess "$RACCESS_BIN" \
    --task-id "$TASK_ID"
