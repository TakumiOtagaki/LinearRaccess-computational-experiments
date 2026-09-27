#!/bin/bash
#$ -cwd
#$ -l gpu_1=1
#$ -l h_rt=06:00:00
#$ -N deepraccess_gpu
#$ -o logs/uge/deepraccess_gpu_$JOB_ID.out
#$ -j y
#
# DeepRaccess (FCN_structured.pth, w=35) via benchmark/scripts/run_deepraccess.py.
# GPU: batch inference on the 12,872 Table 1 sequences and the 20 long ones
# (throughput), then per-sequence (batch size 1) timing on the 52-sequence
# set and the 20 long ones.
set -euo pipefail
PROJECT_ROOT=/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper
cd "${PROJECT_ROOT}"
PY=vendor/DeepRaccess/.venv/bin/python
R=benchmark/scripts/run_deepraccess.py
OUT=benchmark/results/deepraccess_comparison/deepraccess
echo "Start: $(date -u +%FT%TZ) host=$(hostname)"
sha256sum vendor/DeepRaccess/path/FCN_structured.pth > benchmark/results/deepraccess_comparison/checkpoint_sha256.txt
git -C vendor/DeepRaccess rev-parse HEAD > benchmark/results/deepraccess_comparison/deepraccess_commit.txt
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
$PY $R --fasta benchmark/data/bprna_100_1000.fa --ids benchmark/data/bprna_100_1000_list.txt --outdir $OUT/bprna12872 --mode batch --device cuda
$PY $R --fasta benchmark/data/bpRNA_1m_90_filtered_clean.fa --ids benchmark/data/bprna_long_1000_4381_list.txt --outdir $OUT/long20 --mode batch --device cuda
$PY $R --fasta benchmark/data/bprna_100_1000.fa --ids benchmark/data/bprna_52_list.txt --outdir $OUT/timing_gpu_bprna52 --mode single --device cuda
$PY $R --fasta benchmark/data/bpRNA_1m_90_filtered_clean.fa --ids benchmark/data/bprna_long_1000_4381_list.txt --outdir $OUT/timing_gpu_long20 --mode single --device cuda
echo "End: $(date -u +%FT%TZ)"
