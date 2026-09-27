#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_4=1
#$ -l h_rt=03:00:00
#$ -N linsamp_acc_vs_k_r
#$ -j y
#$ -o logs/uge/linsamp_acc_vs_k_reduced_$JOB_ID.out

set -euo pipefail

PROJECT_ROOT="/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper"
cd "${PROJECT_ROOT}"
mkdir -p logs/uge benchmark/results

echo "Host: $(hostname)"
echo "Start: $(date -u +%FT%TZ)"

python3 benchmark/scripts/run_linearsampling_accuracy_vs_k.py \
    --n-per-bin 3 \
    --k-values 50000,1000000,5000000 \
    --out benchmark/results/linearsampling_accuracy_vs_k_reduced.json

echo "End: $(date -u +%FT%TZ)"
