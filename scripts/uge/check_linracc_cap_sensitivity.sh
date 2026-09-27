#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=03:00:00
#$ -N cap_sensitivity
#$ -j y
#$ -o logs/uge/cap_sensitivity_$JOB_ID.out

set -euo pipefail

PROJECT_ROOT="/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper"
cd "${PROJECT_ROOT}"
mkdir -p logs/uge benchmark/results

echo "Host: $(hostname)"
echo "Start: $(date -u +%FT%TZ)"

python3 benchmark/scripts/check_linracc_cap_sensitivity.py

echo "End: $(date -u +%FT%TZ)"
