#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=02:00:00
#$ -N calibrate_lsamp_k
#$ -j y
#$ -o logs/uge/calibrate_lsamp_k_$JOB_ID.out

set -euo pipefail

echo "DEBUG whoami=$(whoami) pwd0=$(pwd) BASH_SOURCE0=${BASH_SOURCE[0]}"
PROJECT_ROOT="/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper"
cd "${PROJECT_ROOT}"
echo "DEBUG pwd1=$(pwd)"

mkdir -p logs/uge benchmark/results

echo "Host: $(hostname)"
echo "Start: $(date -u +%FT%TZ)"

/usr/bin/time -v python3 benchmark/scripts/calibrate_linearsampling_k.py \
    --out benchmark/results/linearsampling_k_calibration.json

echo "End: $(date -u +%FT%TZ)"
