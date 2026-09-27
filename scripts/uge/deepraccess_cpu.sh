#!/bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=06:00:00
#$ -N deepraccess_cpu
#$ -o logs/uge/deepraccess_cpu_$JOB_ID.out
#$ -j y
#
# DeepRaccess (FCN_structured.pth, w=35) via benchmark/scripts/run_deepraccess.py.
# CPU, one thread: per-sequence (batch size 1) timing on the 52-sequence
# set and the 20 long ones, for comparison with the single-core runs of the
# other tools.
set -euo pipefail
PROJECT_ROOT=/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper
cd "${PROJECT_ROOT}"
PY=vendor/DeepRaccess/.venv/bin/python
R=benchmark/scripts/run_deepraccess.py
OUT=benchmark/results/deepraccess_comparison/deepraccess
echo "Start: $(date -u +%FT%TZ) host=$(hostname)"
sha256sum vendor/DeepRaccess/path/FCN_structured.pth > benchmark/results/deepraccess_comparison/checkpoint_sha256.txt
git -C vendor/DeepRaccess rev-parse HEAD > benchmark/results/deepraccess_comparison/deepraccess_commit.txt
mkdir -p $OUT
/usr/bin/time -v -o $OUT/timing_cpu1_bprna52.time $PY $R --fasta benchmark/data/bprna_100_1000.fa --ids benchmark/data/bprna_52_list.txt --outdir $OUT/timing_cpu1_bprna52 --mode single --device cpu --threads 1
/usr/bin/time -v -o $OUT/timing_cpu1_long20.time $PY $R --fasta benchmark/data/bpRNA_1m_90_filtered_clean.fa --ids benchmark/data/bprna_long_1000_4381_list.txt --outdir $OUT/timing_cpu1_long20 --mode single --device cpu --threads 1
echo "End: $(date -u +%FT%TZ)"
