#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=24:00:00
#$ -N bprna_w_cap_sweep

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${PROJECT_ROOT}"

# Item 11: does LinearRaccess's beam/cap approximation still hold at window
# lengths (20, 30, 50 nt) approaching or exceeding the default 30-nt loop
# caps -- the paper's own motivating applications (siRNA/ASO ~15-25nt, RBP
# footprints ~50nt)? Raises c_multi/c_hairpin to 60 (>= the largest tested
# w=50) so capped contexts can in principle still produce a segment event
# long enough to cover these windows; keeps w=1,3,7,10 alongside for
# continuity with Table 1. Runs fresh Raccess Exact too (not reused from
# Table 1/S3): those were computed with access_len=1,3,7,10 only, and the
# w=20/30/50 marginals were never extracted.
#
# Deliberately avoids mktemp/$TMPDIR (see bprna_accuracy_comparison.sh for
# why); every scratch path here is an explicit directory under the repo.
UV_BIN="/home/3/ub07263/.local/bin/uv"

BPRNA_FASTA="${PROJECT_ROOT}/benchmark/data/bpRNA_1m_90_filtered_clean.fa"
SEQ_LIST="${SEQ_LIST:-${PROJECT_ROOT}/benchmark/data/bprna_300_1000_list.txt}"
RESULTS_DIR="${PROJECT_ROOT}/benchmark/results/bprna_w_cap_sweep_300_1000"
WORK_ROOT="${RESULTS_DIR}/.work"
RACCESS_BIN="${PROJECT_ROOT}/vendor/raccess/src/raccess/run_raccess"
LINRACC_BIN="${PROJECT_ROOT}/vendor/LinearRaccess/build-raccess/LinRacc"

ACCESS_LENS="1,3,7,10,20,30,50"
BEAM=100
C_MULTI=60
C_HAIRPIN=60
N_PARALLEL="${N_PARALLEL:-16}"

mkdir -p "${RESULTS_DIR}/exact" "${RESULTS_DIR}/linear_b100_c60" "${WORK_ROOT}"

TOTAL_SEQ=$(wc -l < "${SEQ_LIST}")
START_IDX="${START_IDX:-1}"
END_IDX="${END_IDX:-${TOTAL_SEQ}}"

echo "Processing sequences ${START_IDX} to ${END_IDX} (of ${TOTAL_SEQ}), N_PARALLEL=${N_PARALLEL}"

process_sequence() {
    local LINE_NUM=$1
    local SEQ_LINE=$(sed -n "${LINE_NUM}p" "${SEQ_LIST}")
    local SEQ_ID=$(echo "${SEQ_LINE}" | cut -f1)
    local SEQ_LENGTH=$(echo "${SEQ_LINE}" | cut -f2)

    if [ -f "${RESULTS_DIR}/${SEQ_ID}_summary.json" ]; then
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

    local RACCESS_EXACT_OUT="${RESULTS_DIR}/exact/${SEQ_ID}.txt"
    local LINRACC_OUT="${RESULTS_DIR}/linear_b100_c60/${SEQ_ID}.txt"

    local START_TIME=$(date +%s)

    if /usr/bin/time -v "${RACCESS_BIN}" \
        -seqfile="${SEQ_FILE}" \
        -access_len="${ACCESS_LENS}" \
        -max_span="${SEQ_LENGTH}" \
        -outfile="${RACCESS_EXACT_OUT}" \
        2>"${TEMP_DIR}/time_exact.log"; then
        local EXACT_STATUS="SUCCESS"
        local EXACT_EXIT=0
    else
        local EXACT_STATUS="FAILED"
        local EXACT_EXIT=$?
    fi
    local EXACT_MEM=$(grep "Maximum resident set size" "${TEMP_DIR}/time_exact.log" 2>/dev/null | awk '{print $6}' || echo "null")

    if /usr/bin/time -v "${LINRACC_BIN}" \
        -seqfile="${SEQ_FILE}" \
        -outfile="${LINRACC_OUT}" \
        -access_len="${ACCESS_LENS}" \
        -engine=raccess \
        -energy=turner1999 \
        -beam="${BEAM}" \
        -c_multi="${C_MULTI}" \
        -c_hairpin="${C_HAIRPIN}" \
        -probabilities \
        -metadata \
        2>"${TEMP_DIR}/time_linear.log"; then
        local LINEAR_STATUS="SUCCESS"
        local LINEAR_EXIT=0
    else
        local LINEAR_STATUS="FAILED"
        local LINEAR_EXIT=$?
    fi
    local LINEAR_MEM=$(grep "Maximum resident set size" "${TEMP_DIR}/time_linear.log" 2>/dev/null | awk '{print $6}' || echo "null")

    local END_TIME=$(date +%s)
    local ELAPSED=$((END_TIME - START_TIME))

    "${UV_BIN}" run python3 << PYEOF
import json

results = {
    "seq_id": "${SEQ_ID}",
    "seq_length": ${SEQ_LENGTH},
    "elapsed_sec": ${ELAPSED},
    "access_lens": [1, 3, 7, 10, 20, 30, 50],
    "methods": {
        "raccess-exact": {
            "status": "${EXACT_STATUS}",
            "exit_code": ${EXACT_EXIT},
            "memory_kb": ${EXACT_MEM:-null},
            "max_span": ${SEQ_LENGTH}
        },
        "linearraccess-b100-c60": {
            "status": "${LINEAR_STATUS}",
            "exit_code": ${LINEAR_EXIT},
            "memory_kb": ${LINEAR_MEM:-null},
            "beam": ${BEAM},
            "c_multi": ${C_MULTI},
            "c_hairpin": ${C_HAIRPIN}
        }
    }
}

with open("${RESULTS_DIR}/${SEQ_ID}_summary.json", "w") as f:
    json.dump(results, f, indent=2)
PYEOF

    rm -rf "${TEMP_DIR}"
    echo "Completed ${SEQ_ID} (${SEQ_LENGTH} nt): exact=${EXACT_STATUS}, linear=${LINEAR_STATUS} in ${ELAPSED}s"
}

export -f process_sequence
export BPRNA_FASTA SEQ_LIST RESULTS_DIR WORK_ROOT RACCESS_BIN LINRACC_BIN ACCESS_LENS BEAM C_MULTI C_HAIRPIN UV_BIN

CURRENT_JOBS=0
for LINE_NUM in $(seq "${START_IDX}" "${END_IDX}"); do
    while [ "$(jobs -r | wc -l)" -ge "${N_PARALLEL}" ]; do
        sleep 0.5
    done
    process_sequence "${LINE_NUM}" &
    CURRENT_JOBS=$((CURRENT_JOBS + 1))
    if [ $((CURRENT_JOBS % 100)) -eq 0 ]; then
        echo "Launched ${CURRENT_JOBS} jobs so far..."
    fi
done

wait
echo "Completed sequences ${START_IDX}-${END_IDX}"
