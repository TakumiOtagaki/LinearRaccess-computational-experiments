#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=02:00:00
#$ -N verify_k_crossover
#$ -j y
#$ -o logs/uge/verify_k_crossover_$JOB_ID.out

set -euo pipefail

PROJECT_ROOT="/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper"
cd "${PROJECT_ROOT}"
mkdir -p logs/uge benchmark/results

echo "Host: $(hostname)"
echo "Start: $(date -u +%FT%TZ)"

python3 benchmark/scripts/verify_k_crossover.py

echo "End: $(date -u +%FT%TZ)"
