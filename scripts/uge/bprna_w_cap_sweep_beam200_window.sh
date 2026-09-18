#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=24:00:00
#$ -N bprna_w_cap_sweep_b200_win

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${PROJECT_ROOT}"

# Follow-up to bprna_w_cap_sweep.sh (beam=100, c=60): author's hypothesis is
# that beam=100 may be insufficient at large w, since a long unpaired run
# has a lower Boltzmann weight than typical paired states and could be
# pruned from the beam before the aggregation ever sees it -- separate from
# whether the loop cap is large enough. Tests this with beam=200 at the same
# cap=60 and access_lens, reusing the exact/ reference already computed by
# bprna_w_cap_sweep.sh (unaffected by beam/cap, so no need to rerun it).
# Also runs Raccess Window (max_span=100) at the same access_lens for the
# comparison the author asked for -- independent of beam/cap entirely.
UV_BIN="/home/3/ub07263/.local/bin/uv"

BPRNA_FASTA="${PROJECT_ROOT}/benchmark/data/bpRNA_1m_90_filtered_clean.fa"
SEQ_LIST="${SEQ_LIST:-${PROJECT_ROOT}/benchmark/data/bprna_300_1000_list.txt}"
RESULTS_DIR="${PROJECT_ROOT}/benchmark/results/bprna_w_cap_sweep_300_1000"
WORK_ROOT="${RESULTS_DIR}/.work_b200win"
RACCESS_BIN="${PROJECT_ROOT}/vendor/raccess/src/raccess/run_raccess"
LINRACC_BIN="${PROJECT_ROOT}/vendor/LinearRaccess/build-raccess/LinRacc"

ACCESS_LENS="1,3,7,10,20,30,50"
BEAM=200
C_MULTI=60
C_HAIRPIN=60
N_PARALLEL="${N_PARALLEL:-16}"

mkdir -p "${RESULTS_DIR}/linear_b200_c60" "${RESULTS_DIR}/window_w100" "${WORK_ROOT}"

TOTAL_SEQ=$(wc -l < "${SEQ_LIST}")
START_IDX="${START_IDX:-1}"
END_IDX="${END_IDX:-${TOTAL_SEQ}}"

echo "Processing sequences ${START_IDX} to ${END_IDX} (of ${TOTAL_SEQ}), N_PARALLEL=${N_PARALLEL}"

process_sequence() {
    local LINE_NUM=$1
    local SEQ_LINE=$(sed -n "${LINE_NUM}p" "${SEQ_LIST}")
    local SEQ_ID=$(echo "${SEQ_LINE}" | cut -f1)
    local SEQ_LENGTH=$(echo "${SEQ_LINE}" | cut -f2)

    if [ -f "${RESULTS_DIR}/${SEQ_ID}_b200win_summary.json" ]; then
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

    local LINRACC_OUT="${RESULTS_DIR}/linear_b200_c60/${SEQ_ID}.txt"
    local WINDOW_OUT="${RESULTS_DIR}/window_w100/${SEQ_ID}.txt"

    local START_TIME=$(date +%s)

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

    if /usr/bin/time -v "${RACCESS_BIN}" \
        -seqfile="${SEQ_FILE}" \
        -access_len="${ACCESS_LENS}" \
        -max_span=100 \
        -outfile="${WINDOW_OUT}" \
        2>"${TEMP_DIR}/time_window.log"; then
        local WINDOW_STATUS="SUCCESS"
        local WINDOW_EXIT=0
    else
        local WINDOW_STATUS="FAILED"
        local WINDOW_EXIT=$?
    fi

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
        "linearraccess-b200-c60": {"status": "${LINEAR_STATUS}", "exit_code": ${LINEAR_EXIT}, "beam": ${BEAM}, "c_multi": ${C_MULTI}, "c_hairpin": ${C_HAIRPIN}},
        "raccess-w100": {"status": "${WINDOW_STATUS}", "exit_code": ${WINDOW_EXIT}, "max_span": 100}
    }
}
with open("${RESULTS_DIR}/${SEQ_ID}_b200win_summary.json", "w") as f:
    json.dump(results, f, indent=2)
PYEOF

    rm -rf "${TEMP_DIR}"
    echo "Completed ${SEQ_ID} (${SEQ_LENGTH} nt): linear_b200=${LINEAR_STATUS}, window=${WINDOW_STATUS} in ${ELAPSED}s"
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
