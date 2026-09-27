#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=12:00:00
#$ -N linracc_em_beam
#$ -j y
#$ -o logs/uge/linracc_em_beam_$JOB_ID.out

set -euo pipefail
cd /gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper
echo "Host: $(hostname)  Start: $(date -u +%FT%TZ)"
python3 benchmark/scripts/run_linracc_energy_model_beam_check.py
echo "End: $(date -u +%FT%TZ)"
