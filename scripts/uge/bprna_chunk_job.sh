#!/usr/bin/env bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=04:00:00
#$ -N linearraccess_bprna_chunk
#$ -t 1-27
#$ -o logs/chunk_$JOB_ID.$TASK_ID.out
#$ -e logs/chunk_$JOB_ID.$TASK_ID.err
#$ -tc 10

set -euo pipefail

REPO_ROOT="${BPRNA_REPO_ROOT:-/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper}"
TASK_ID="${SGE_TASK_ID:?SGE_TASK_ID is required}"
CHUNK_DIR="${BPRNA_CHUNK_DIR:-$REPO_ROOT/benchmark/data/bprna_chunks}"
OUTDIR="${BPRNA_OUTDIR:-$REPO_ROOT/benchmark/results/bprna_runs/${JOB_ID:-manual}}"
LINRACC_BIN="${BPRNA_LINRACC:-$REPO_ROOT/vendor/LinearRaccess/build-raccess/LinRacc}"
CHUNK_FILE="$CHUNK_DIR/chunk_$(printf '%02d' "$TASK_ID").fa"

ARGS=(--chunk "$CHUNK_FILE" --outdir "$OUTDIR" --linracc "$LINRACC_BIN"
      --task-id "$TASK_ID" --parallel "${BPRNA_PARALLEL:-10}"
      --timeout "${BPRNA_TIMEOUT:-3600}")
if [[ "${BPRNA_NO_FAST_LOGSUMEXP:-0}" == 1 ]]; then
    ARGS+=(--no-fast-logsumexp)
fi
exec python3 "$REPO_ROOT/benchmark/scripts/run_bprna_chunk.py" "${ARGS[@]}"
