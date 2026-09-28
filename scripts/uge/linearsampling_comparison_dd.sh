#!/bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=24:00:00
#$ -N ls_cmp_dd
#$ -o logs/uge/linearsampling_comparison_dd_$JOB_ID.out
#$ -j y
#
# Rerun of the LinearRaccess side (exact reference + LinearRaccess configs) of
# the LinearSampling comparison with the double-double build; LinearSampling
# caches are copied in beforehand and reused.  Env: SET (bprna52|long),
# PHASE (accuracy|timing), WORKDIR (qsub -v).
set -euo pipefail
PROJECT_ROOT=/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper
cd "${PROJECT_ROOT}"
: "${SET:?}" "${PHASE:?}" "${WORKDIR:?}"
JOBS=16; [ "${SET}" = long ] && JOBS=10
echo "Start: $(date -u +%FT%TZ) host=$(hostname) set=${SET} phase=${PHASE}"
python3 benchmark/scripts/run_linearsampling_comparison_capb.py --set "${SET}" --phase "${PHASE}" \
    --work-dir "${WORKDIR}" --jobs "${JOBS}"
echo "End: $(date -u +%FT%TZ)"
