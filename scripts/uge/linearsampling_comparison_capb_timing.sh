#!/bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=12:00:00
#$ -N ls_cmp_timing
#$ -o logs/uge/linearsampling_comparison_capb_timing_$JOB_ID.out
#$ -j y
#
# LinearRaccess vs LinearSampling, timing phase: see the docstring of
# benchmark/scripts/run_linearsampling_comparison_capb.py.  Uses the frozen
# LinRacc in benchmark/results/linearsampling_comparison_capb/bin/.
set -euo pipefail
PROJECT_ROOT=/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper
cd "${PROJECT_ROOT}"
echo "Start: $(date -u +%FT%TZ) host=$(hostname)"
python3 benchmark/scripts/run_linearsampling_comparison_capb.py --phase timing
echo "End: $(date -u +%FT%TZ)"
