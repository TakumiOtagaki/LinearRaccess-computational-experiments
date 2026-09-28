#!/bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=24:00:00
#$ -N byloop_cmp
#$ -o logs/uge/byloop_comparison_$JOB_ID.out
#$ -j y
#
# Loop-type-resolved accessibility, LinearRaccess vs LinearSampling
# (benchmark/scripts/run_byloop_comparison.py).  Env: SET (bprna52|long),
# PHASE (accuracy|timing) (qsub -v).
set -euo pipefail
PROJECT_ROOT=/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper
cd "${PROJECT_ROOT}"
: "${SET:?}" "${PHASE:?}"
JOBS=16; [ "${SET}" = long ] && JOBS=10
echo "Start: $(date -u +%FT%TZ) host=$(hostname) set=${SET} phase=${PHASE}"
python3 benchmark/scripts/run_byloop_comparison.py --set "${SET}" --phase "${PHASE}" \
    --work-dir "benchmark/results/byloop_comparison_${SET}" --jobs "${JOBS}"
echo "End: $(date -u +%FT%TZ)"
