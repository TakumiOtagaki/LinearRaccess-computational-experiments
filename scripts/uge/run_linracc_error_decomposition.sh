#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=12:00:00
#$ -N linracc_err_decomp
#$ -j y
#$ -o logs/uge/linracc_err_decomp_$JOB_ID.out

set -euo pipefail
PROJECT_ROOT="/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper"
cd "${PROJECT_ROOT}"
mkdir -p logs/uge
echo "Host: $(hostname)  Start: $(date -u +%FT%TZ)"
python3 benchmark/scripts/run_linracc_error_decomposition.py --jobs 16
echo "End: $(date -u +%FT%TZ)"
