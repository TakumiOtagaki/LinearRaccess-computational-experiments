#!/bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=24:00:00
#$ -N ls_long_timing
#$ -o logs/uge/linearsampling_comparison_long_timing_$JOB_ID.out
#$ -j y
#
# LinearRaccess vs LinearSampling on 20 bpRNA sequences of 1,000-3,541 nt
# (longer than the beam, so pruning is active), timing phase; see
# run_linearsampling_comparison_capb.py --set long.  Frozen LinRacc (5b0c139)
# in benchmark/results/linearsampling_comparison_long/bin/.
set -euo pipefail
PROJECT_ROOT=/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper
cd "${PROJECT_ROOT}"
echo "Start: $(date -u +%FT%TZ) host=$(hostname)"
python3 benchmark/scripts/run_linearsampling_comparison_capb.py --set long --phase timing
echo "End: $(date -u +%FT%TZ)"
