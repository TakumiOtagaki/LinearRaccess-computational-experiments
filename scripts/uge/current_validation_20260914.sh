#!/usr/bin/env bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=04:00:00
#$ -N current_validation
#$ -t 1-2
#$ -o benchmark/scripts/uge/logs/current_validation_$JOB_ID.$TASK_ID.out
#$ -e benchmark/scripts/uge/logs/current_validation_$JOB_ID.$TASK_ID.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ROOT=/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper
cd "$ROOT"
OUT=benchmark/results/current_validation_20260914
mkdir -p "$OUT/runtime/$JOB_ID.$SGE_TASK_ID"
BIN="$ROOT/$OUT/runtime/$JOB_ID.$SGE_TASK_ID/LinRacc"
cp vendor/LinearRaccess/build-raccess/LinRacc "$BIN"
if [[ "$SGE_TASK_ID" == 2 ]]; then
    exec timeout 3h python3 benchmark/scripts/run_scaling.py \
        --manifest benchmark/manifests/e3_current_20260914.json --tier long \
        --software-dir vendor/LinearRaccess --linracc "$BIN" \
        --output-dir "$OUT/e3"
fi
RBIN="$ROOT/$OUT/runtime/$JOB_ID.$SGE_TASK_ID/run_raccess"
cp vendor/Raccess/src/raccess/run_raccess "$RBIN"
for sample in regression bprna; do
    if [[ "$sample" == regression ]]; then
        INPUT=benchmark/data/e1_raccess_inputs.fa
        MANIFEST=benchmark/manifests/e2_approximation_envelope.json
    else
        INPUT=benchmark/data/bprna_reference_20260914.fa
        MANIFEST=benchmark/manifests/e2_bprna_reference_20260914.json
    fi
    timeout 90m python3 benchmark/scripts/run_raccess_compare.py \
        --bundle "$OUT/e1_$sample" --input "$INPUT" \
        --raccess-binary "$RBIN" --linracc-binary "$BIN" \
        --raccess-repo vendor/Raccess \
        --raccess-commit ad29d098b62606aa4051951469780b4e7be54536 \
        --raccess-patch-note 'Existing binary frozen by SHA256; local checkout has def_opt_main.hpp changes and an added energy_model_api.hpp. Build-to-source correspondence is not independently reconstructed.'
    timeout 90m python3 benchmark/scripts/run_algorithm_validation.py \
        --manifest "$MANIFEST" --software-dir vendor/LinearRaccess \
        --linracc "$BIN" --output-dir "$OUT/e2_$sample"
done
