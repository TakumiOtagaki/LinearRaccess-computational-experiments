#!/usr/bin/env bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=04:00:00
#$ -N bprna_memory
#$ -t 1-27
#$ -o logs/bprna_memory_$JOB_ID.$TASK_ID.out
#$ -e logs/bprna_memory_$JOB_ID.$TASK_ID.err
#$ -tc 10

set -euo pipefail

REPO_ROOT="/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper"
TASK_ID="${SGE_TASK_ID:?SGE_TASK_ID is required}"
OUTDIR="${REPO_ROOT}/benchmark/results/bprna_memory_runs/${JOB_ID:-manual}"
CHUNK_DIR="${REPO_ROOT}/benchmark/data/bprna_chunks"
LINRACC_BIN="${REPO_ROOT}/vendor/LinearRaccess/build-raccess/LinRacc"
CHUNK_FILE="$CHUNK_DIR/chunk_$(printf '%02d' "$TASK_ID").fa"

mkdir -p "$OUTDIR"

cd "$REPO_ROOT"
exec /home/3/ub07263/.local/bin/uv run python "$REPO_ROOT/benchmark/scripts/run_bprna_chunk_with_memory.py" \
    --chunk "$CHUNK_FILE" \
    --outdir "$OUTDIR" \
    --linracc "$LINRACC_BIN" \
    --task-id "$TASK_ID" \
    --parallel 10 \
    --timeout 3600
