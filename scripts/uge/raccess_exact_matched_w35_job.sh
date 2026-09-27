#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=12:00:00
#$ -t 1-20:1
#$ -N racc_exact_w35
#$ -o logs/raccess_exact_matched_$JOB_ID.$TASK_ID.out
#$ -e logs/raccess_exact_matched_$JOB_ID.$TASK_ID.err

set -euo pipefail

REPO_ROOT="/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper"
TASK_ID="${SGE_TASK_ID:?SGE_TASK_ID is required}"
BATCH_NUM=$(printf "%02d" "$TASK_ID")
BATCH_FILE="${REPO_ROOT}/benchmark/data/raccess_exact_batches/batch_${BATCH_NUM}.txt"
FASTA="${REPO_ROOT}/benchmark/data/bpRNA_1m_90_filtered_clean.fa"
RACCESS_BIN="${REPO_ROOT}/vendor/raccess/src/raccess/run_raccess"
OUTDIR="${REPO_ROOT}/benchmark/results/raccess_exact_matched_w35/${JOB_ID:-manual}"

mkdir -p "$OUTDIR" "${REPO_ROOT}/logs"

cd "$REPO_ROOT"
exec python3 "${REPO_ROOT}/benchmark/scripts/run_raccess_exact_chunk_matched.py" \
    --batch "$BATCH_FILE" \
    --fasta "$FASTA" \
    --outdir "$OUTDIR" \
    --raccess "$RACCESS_BIN" \
    --task-id "$TASK_ID" \
    --access-len 1,3,7,10,35
