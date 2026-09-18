#!/bin/bash
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l cpu_16=1
#$ -N e3_scaling_linearraccess
#$ -o benchmark/scripts/uge/logs/e3_scaling_$JOB_ID.$TASK_ID.out
#$ -e benchmark/scripts/uge/logs/e3_scaling_$JOB_ID.$TASK_ID.err
#$ -j y

set -euo pipefail

REPO_ROOT="/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper"
cd $REPO_ROOT

echo "=== E3 Scaling Benchmark ==="
echo "Date: $(date)"
echo "Beam: 100 (updated from 200)"
echo "Prefix lengths: 1000, 2000, 4000, 8000, 16000, 29903"
echo "Repetitions: 3"
echo ""

python benchmark/scripts/run_scaling.py \
  --manifest benchmark/manifests/scaling_sars2_long.draft.json \
  --tier long \
  --software-dir vendor/LinearRaccess \
  --linracc vendor/LinearRaccess/build-raccess/LinRacc \
  --output-dir benchmark/results/e3_long_scaling_2026-09-10

echo "=== Done ==="
echo "Exit code: $?"
