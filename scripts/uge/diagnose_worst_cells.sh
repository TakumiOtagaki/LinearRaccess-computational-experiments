#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=06:00:00
#$ -N diag_worst_cells
#$ -j y
#$ -o logs/uge/diag_worst_cells_$JOB_ID.out

set -euo pipefail
PROJECT_ROOT="/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper"
cd "${PROJECT_ROOT}"
mkdir -p logs/uge benchmark/results
echo "Host: $(hostname)  Start: $(date -u +%FT%TZ)"
python3 benchmark/scripts/diagnose_worst_cells.py
echo "End: $(date -u +%FT%TZ)"
