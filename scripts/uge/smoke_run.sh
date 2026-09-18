#!/bin/bash
# UGE (SGE) compute-node smoke run for LinearRaccess paper scaling runners.
# Site: xxx-cluster (UGE); paths are absolute. Not TSUBAME pjsub.
# Run both synthetic smoke tiers and print short provenance for the log.

#$ -N linracc-smoke
#$ -l h_rt=1:00:00
#$ -l m_mem_free=4G
#$ -cwd

set -euo pipefail

PAPER=/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper
SDIR="$PAPER/vendor/LinearRaccess"
LINRACC="$SDIR/build/LinRacc"
SMOKE_OUT="$PAPER/scratch/smoke"

echo "HOSTNAME=$(hostname)"
echo "NPROC=$(nproc)"
echo "PAPER_COMMIT=$(git -C "$PAPER" rev-parse HEAD)"
echo "SOFTWARE_COMMIT=$(git -C "$SDIR" rev-parse HEAD)"
echo "GCC=$(gcc --version | head -1)"
echo "PYTHON=$(python3 --version)"

rm -rf "$SMOKE_OUT"
mkdir -p "$SMOKE_OUT"

python3 "$PAPER/benchmark/scripts/run_scaling.py" \
  --manifest "$PAPER/benchmark/manifests/scaling_smoke.json" \
  --tier smoke \
  --software-dir "$SDIR" \
  --linracc "$LINRACC" \
  --output-dir "$SMOKE_OUT/scaling-smoke-uge-v1"

python3 "$PAPER/benchmark/scripts/run_work_scaling.py" \
  --manifest "$PAPER/benchmark/manifests/scaling_smoke.json" \
  --tier smoke \
  --software-dir "$SDIR" \
  --linracc "$LINRACC" \
  --output-dir "$SMOKE_OUT/work-scaling-smoke-uge-v1"

echo "SMOKE_OK"