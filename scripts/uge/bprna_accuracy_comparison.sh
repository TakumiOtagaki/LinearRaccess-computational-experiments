#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=24:00:00
#$ -t 1-27:1
#$ -tc 27
#$ -N bprna_acc_cmp
#$ -o logs/uge/bprna_acc_$TASK_ID.out
#$ -e logs/uge/bprna_acc_$TASK_ID.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${PROJECT_ROOT}"

# Deliberately never touch $TMPDIR/mktemp: on this cluster, something outside
# this script (not this repo, not the raccess/LinRacc binaries) tries to
# mkdir a scheduler-side scratch path (seen as /var/spool/benchmark and
# /var/spool/tmp) and aborts the job before it even starts, for reasons we
# could not isolate. Every scratch path here is an explicit directory under
# the repo instead, so nothing in this script depends on TMPDIR at all.
UV_BIN="/home/3/ub07263/.local/bin/uv"

# Configuration
BPRNA_FASTA="${PROJECT_ROOT}/benchmark/data/bpRNA_1m_90_filtered_clean.fa"
SEQ_LIST="${PROJECT_ROOT}/benchmark/data/bprna_100_1000_list.txt"
RESULTS_DIR="${PROJECT_ROOT}/benchmark/results/bprna_accuracy_comparison"
WORK_ROOT="${RESULTS_DIR}/.work"
LOGDIR="${PROJECT_ROOT}/logs/uge"
RACCESS_BIN="${PROJECT_ROOT}/vendor/raccess/src/raccess/run_raccess"
LINRACC_BIN="${PROJECT_ROOT}/vendor/LinearRaccess/build-raccess/LinRacc"

# Settings
ACCESS_LENS="1,3,7,10"
BEAM=100
C_MULTI=30
C_HAIRPIN=30
N_PARALLEL="${N_PARALLEL:-10}"

# Create output directories
mkdir -p "${RESULTS_DIR}/exact"
mkdir -p "${RESULTS_DIR}/linear_b100"
mkdir -p "${RESULTS_DIR}/window_w100"
mkdir -p "${WORK_ROOT}"
mkdir -p "${LOGDIR}"

# Calculate sequence range for this task
TOTAL_SEQ=12872
if [ -n "${SGE_TASK_ID:-}" ] && [ "${SGE_TASK_ID}" != "undefined" ]; then
    NUM_TASKS=27
    BATCH_SIZE=$(( (TOTAL_SEQ + NUM_TASKS - 1) / NUM_TASKS ))  # Ceiling division
    START_IDX=$(( (SGE_TASK_ID - 1) * BATCH_SIZE + 1 ))
    END_IDX=$(( SGE_TASK_ID * BATCH_SIZE ))
    if [ ${END_IDX} -gt ${TOTAL_SEQ} ]; then
        END_IDX=${TOTAL_SEQ}
    fi
else
    # No SGE array: process the full list in one invocation. Combined with
    # the already-done skip in process_sequence, this run can be resumed by
    # re-invoking the script (e.g. across separate interactive allocations).
    START_IDX="${START_IDX:-1}"
    END_IDX="${END_IDX:-${TOTAL_SEQ}}"
fi

echo "Processing sequences ${START_IDX} to ${END_IDX} (of ${TOTAL_SEQ}), N_PARALLEL=${N_PARALLEL}"

# Function to process a single sequence
process_sequence() {
    local LINE_NUM=$1
    local SEQ_LINE=$(sed -n "${LINE_NUM}p" "${SEQ_LIST}")
    local SEQ_ID=$(echo "${SEQ_LINE}" | cut -f1)
    local SEQ_LENGTH=$(echo "${SEQ_LINE}" | cut -f2)

    if [ -f "${RESULTS_DIR}/${SEQ_ID}_summary.json" ]; then
        echo "Skipping ${SEQ_ID} (already done)"
        return 0
    fi

    local TEMP_DIR="${WORK_ROOT}/${SEQ_ID}"
    mkdir -p "${TEMP_DIR}"
    local SEQ_FILE="${TEMP_DIR}/seq.fa"
    
    # Extract sequence
    "${UV_BIN}" run python3 << EOF
from Bio import SeqIO
for record in SeqIO.parse("${BPRNA_FASTA}", "fasta"):
    if record.id == "${SEQ_ID}":
        SeqIO.write([record], "${SEQ_FILE}", "fasta")
        break
EOF
    
    # Output files
    local RACCESS_EXACT_OUT="${RESULTS_DIR}/exact/${SEQ_ID}.txt"
    local RACCESS_EXACT_JSON="${RESULTS_DIR}/exact/${SEQ_ID}.json"
    local LINRACC_OUT="${RESULTS_DIR}/linear_b100/${SEQ_ID}.txt"
    local LINRACC_JSON="${RESULTS_DIR}/linear_b100/${SEQ_ID}.json"
    local RACCESS_W100_OUT="${RESULTS_DIR}/window_w100/${SEQ_ID}.txt"
    local RACCESS_W100_JSON="${RESULTS_DIR}/window_w100/${SEQ_ID}.json"
    
    local START_TIME=$(date +%s)
    
    # 1. Run Raccess exact
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
    
    # 2. Run LinearRaccess
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
    
    # 3. Run Raccess Window
    if /usr/bin/time -v "${RACCESS_BIN}" \
        -seqfile="${SEQ_FILE}" \
        -access_len="${ACCESS_LENS}" \
        -max_span=100 \
        -outfile="${RACCESS_W100_OUT}" \
        2>"${TEMP_DIR}/time_window.log"; then
        local WINDOW_STATUS="SUCCESS"
        local WINDOW_EXIT=0
    else
        local WINDOW_STATUS="FAILED"
        local WINDOW_EXIT=$?
    fi
    local WINDOW_MEM=$(grep "Maximum resident set size" "${TEMP_DIR}/time_window.log" 2>/dev/null | awk '{print $6}' || echo "null")
    
    local END_TIME=$(date +%s)
    local ELAPSED=$((END_TIME - START_TIME))
    
    # Save JSON results
    "${UV_BIN}" run python3 << PYEOF
import json

results = {
    "seq_id": "${SEQ_ID}",
    "seq_length": ${SEQ_LENGTH},
    "elapsed_sec": ${ELAPSED},
    "methods": {
        "raccess-exact": {
            "status": "${EXACT_STATUS}",
            "exit_code": ${EXACT_EXIT},
            "memory_kb": ${EXACT_MEM:-null},
            "max_span": ${SEQ_LENGTH}
        },
        "linearraccess-b100": {
            "status": "${LINEAR_STATUS}",
            "exit_code": ${LINEAR_EXIT},
            "memory_kb": ${LINEAR_MEM:-null},
            "beam": ${BEAM}
        },
        "raccess-w100": {
            "status": "${WINDOW_STATUS}",
            "exit_code": ${WINDOW_EXIT},
            "memory_kb": ${WINDOW_MEM:-null},
            "max_span": 100
        }
    }
}

with open("${RESULTS_DIR}/${SEQ_ID}_summary.json", "w") as f:
    json.dump(results, f, indent=2)
PYEOF
    
    # Cleanup
    rm -rf "${TEMP_DIR}"
    
    echo "Completed ${SEQ_ID} (${SEQ_LENGTH} nt): exact=${EXACT_STATUS}, linear=${LINEAR_STATUS}, window=${WINDOW_STATUS} in ${ELAPSED}s"
}

export -f process_sequence
export BPRNA_FASTA SEQ_LIST RESULTS_DIR WORK_ROOT RACCESS_BIN LINRACC_BIN ACCESS_LENS BEAM C_MULTI C_HAIRPIN UV_BIN

# Run sequences in parallel (max 10 concurrent)
echo "Starting parallel processing with max ${N_PARALLEL} concurrent jobs..."

CURRENT_JOBS=0
for LINE_NUM in $(seq ${START_IDX} ${END_IDX}); do
    # Wait if we have too many parallel jobs
    while [ $(jobs -r | wc -l) -ge ${N_PARALLEL} ]; do
        sleep 0.5
    done
    
    process_sequence ${LINE_NUM} &
    CURRENT_JOBS=$((CURRENT_JOBS + 1))
    
    # Progress report every 10 sequences
    if [ $((CURRENT_JOBS % 10)) -eq 0 ]; then
        echo "Launched ${CURRENT_JOBS} jobs so far..."
    fi
done

# Wait for all background jobs to complete
echo "Waiting for all jobs to complete..."
wait

echo "Completed sequences ${START_IDX}-${END_IDX}"
exit 0
