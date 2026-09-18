#!/usr/bin/env bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=01:00:00
#$ -N bprna_lse_tradeoff
#$ -o benchmark/scripts/uge/logs/lse_tradeoff_$JOB_ID.out
#$ -e benchmark/scripts/uge/logs/lse_tradeoff_$JOB_ID.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 benchmark/scripts/run_bprna_logsumexp_tradeoff.py run \
    --out "${BPRNA_TRADEOFF_OUTDIR:-benchmark/results/bprna_logsumexp_tradeoff_20260914}"
