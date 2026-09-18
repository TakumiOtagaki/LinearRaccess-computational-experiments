#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=24:00:00
#$ -N bprna_beam_sweep

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${PROJECT_ROOT}"

# LinearRaccess-only beam-width sweep for a fixed sequence subset: how does
# accuracy against the exact reference change as beam width increases?
# Reuses the exact/ reference already computed in
# benchmark/results/bprna_accuracy_comparison/exact/ -- this script only
# runs LinearRaccess at BEAM, it does not recompute exact or Raccess Window.
#
# Deliberately avoids mktemp/$TMPDIR (see bprna_accuracy_comparison.sh for
# why); every scratch path here is an explicit directory under the repo.
UV_BIN="/home/3/ub07263/.local/bin/uv"

BPRNA_FASTA="${PROJECT_ROOT}/benchmark/data/bpRNA_1m_90_filtered_clean.fa"
SEQ_LIST="${SEQ_LIST:-${PROJECT_ROOT}/benchmark/data/bprna_300_1000_list.txt}"
BEAM="${BEAM:?set BEAM=50|200|300|...}"
RESULTS_DIR="${PROJECT_ROOT}/benchmark/results/bprna_beam_sweep_300_1000/linear_b${BEAM}"
WORK_ROOT="${RESULTS_DIR}/.work"
LINRACC_BIN="${PROJECT_ROOT}/vendor/LinearRaccess/build-raccess/LinRacc"

ACCESS_LENS="1,3,7,10"
C_MULTI=30
C_HAIRPIN=30
N_PARALLEL="${N_PARALLEL:-16}"

mkdir -p "${RESULTS_DIR}" "${WORK_ROOT}"

TOTAL_SEQ=$(wc -l < "${SEQ_LIST}")
START_IDX="${START_IDX:-1}"
END_IDX="${END_IDX:-${TOTAL_SEQ}}"

echo "beam=${BEAM}: processing sequences ${START_IDX} to ${END_IDX} (of ${TOTAL_SEQ}), N_PARALLEL=${N_PARALLEL}"

process_sequence() {
    local LINE_NUM=$1
    local SEQ_LINE=$(sed -n "${LINE_NUM}p" "${SEQ_LIST}")
    local SEQ_ID=$(echo "${SEQ_LINE}" | cut -f1)
    local SEQ_LENGTH=$(echo "${SEQ_LINE}" | cut -f2)

    if [ -f "${RESULTS_DIR}/${SEQ_ID}.txt" ]; then
        return 0
    fi

    local TEMP_DIR="${WORK_ROOT}/${SEQ_ID}"
    mkdir -p "${TEMP_DIR}"
    local SEQ_FILE="${TEMP_DIR}/seq.fa"

    "${UV_BIN}" run python3 << EOF
from Bio import SeqIO
for record in SeqIO.parse("${BPRNA_FASTA}", "fasta"):
    if record.id == "${SEQ_ID}":
        SeqIO.write([record], "${SEQ_FILE}", "fasta")
        break
EOF

    "${LINRACC_BIN}" \
        -seqfile="${SEQ_FILE}" \
        -outfile="${RESULTS_DIR}/${SEQ_ID}.txt" \
        -access_len="${ACCESS_LENS}" \
        -engine=raccess \
        -energy=turner1999 \
        -beam="${BEAM}" \
        -c_multi="${C_MULTI}" \
        -c_hairpin="${C_HAIRPIN}" \
        -probabilities \
        -metadata

    rm -rf "${TEMP_DIR}"
    echo "done ${SEQ_ID} (${SEQ_LENGTH} nt) beam=${BEAM}"
}

export -f process_sequence
export BPRNA_FASTA SEQ_LIST RESULTS_DIR WORK_ROOT LINRACC_BIN ACCESS_LENS BEAM C_MULTI C_HAIRPIN UV_BIN

CURRENT_JOBS=0
for LINE_NUM in $(seq "${START_IDX}" "${END_IDX}"); do
    while [ "$(jobs -r | wc -l)" -ge "${N_PARALLEL}" ]; do
        sleep 0.5
    done
    process_sequence "${LINE_NUM}" &
    CURRENT_JOBS=$((CURRENT_JOBS + 1))
    if [ $((CURRENT_JOBS % 100)) -eq 0 ]; then
        echo "launched ${CURRENT_JOBS} jobs so far..."
    fi
done

wait
echo "beam=${BEAM} completed sequences ${START_IDX}-${END_IDX}"
