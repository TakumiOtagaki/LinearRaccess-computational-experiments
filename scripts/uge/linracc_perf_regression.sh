#!/bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=3:00:00
#$ -j y
#$ -o logs/uge/linracc_perf_regression_$JOB_ID.out
# Byte-identical output check: LinRacc 59b692c (reference) vs the
# single-pass-accessibility candidate, over a grid of CLI settings.
set -euo pipefail
PROJECT_ROOT=/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper
W=$PROJECT_ROOT/benchmark/results/.work/linracc_perf_regression
cd "$PROJECT_ROOT"
echo "Start: $(date -u +%FT%TZ) host=$(hostname)"
python3 vendor/LinearRaccess/tests/compare_cli_outputs.py \
  --reference "$W/LinRacc_baseline_59b692c" --candidate "$W/LinRacc_candidate" \
  --fasta "$W/reg.fa" --with-raccess --jobs 16
echo "End: $(date -u +%FT%TZ)"
