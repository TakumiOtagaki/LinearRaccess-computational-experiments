#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=24:00:00
#$ -N beam_sweep_exact
#$ -j y
#$ -o logs/uge/beam_sweep_exact_$JOB_ID.out

set -euo pipefail
cd /gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper
echo "Host: $(hostname)  Start: $(date -u +%FT%TZ)"
python3 benchmark/scripts/run_beam_sweep_vs_exact.py
echo "End: $(date -u +%FT%TZ)"
