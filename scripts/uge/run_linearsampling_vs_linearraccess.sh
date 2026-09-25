#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=04:00:00
#$ -N linsamp_vs_linracc
#$ -j y
#$ -o logs/uge/linsamp_vs_linracc_$JOB_ID.out

set -euo pipefail

PROJECT_ROOT="/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper"
cd "${PROJECT_ROOT}"
mkdir -p logs/uge benchmark/results

echo "Host: $(hostname)"
echo "Start: $(date -u +%FT%TZ)"

/usr/bin/time -v python3 benchmark/scripts/run_linearsampling_vs_linearraccess.py \
    --n-per-bin 13 \
    --out benchmark/results/linearsampling_accuracy.json

echo "End: $(date -u +%FT%TZ)"
