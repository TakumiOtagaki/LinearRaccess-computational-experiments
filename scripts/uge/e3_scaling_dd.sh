#!/bin/bash
#$ -cwd
#$ -l h_rt=12:00:00
#$ -l cpu_16=1
#$ -N e3_scaling_dd
#$ -o logs/uge/e3_scaling_dd_$JOB_ID.out
#$ -j y
#
# Section 3.4 SARS-CoV-2 runtime/memory (and its prefix-length scaling),
# rerun for the single-pass accessibility build with loop caps tied to the
# beam and w in {1,3,7,10,35}.  Same protocol as e3_scaling_run.sh:
# run_scaling.py, beam=200 (the new default), Turner 2004, one warm-up + three measured
# repetitions per prefix length, one process at a time on a 16-core slot.
# run_scaling.py requires a clean paper repository and the LinearRaccess
# checkout at software.lock.json's commit.

set -euo pipefail
REPO_ROOT=/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper
cd "${REPO_ROOT}"
echo "Start: $(date -u +%FT%TZ) host=$(hostname)"
python3 benchmark/scripts/run_scaling.py \
  --manifest benchmark/manifests/scaling_sars2_long_b200.json \
  --tier long \
  --software-dir vendor/LinearRaccess \
  --linracc vendor/LinearRaccess/build-raccess/LinRacc \
  --output-dir benchmark/results/e3_long_scaling_dd_2026-09-28
echo "End: $(date -u +%FT%TZ)"
