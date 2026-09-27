#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=12:00:00
#$ -N exact_ref_si
#$ -j y
#$ -o logs/uge/exact_ref_si_$JOB_ID.out

set -euo pipefail
PROJECT_ROOT="/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper"
cd "${PROJECT_ROOT}"
mkdir -p logs/uge benchmark/results
echo "Host: $(hostname)  Start: $(date -u +%FT%TZ)"
python3 benchmark/scripts/run_exact_reference_si.py --n-per-bin 13 --jobs 16
echo "End: $(date -u +%FT%TZ)"
