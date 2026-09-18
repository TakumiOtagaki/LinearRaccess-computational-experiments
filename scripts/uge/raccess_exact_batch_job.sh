#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=12:00:00
#$ -t 1-20:1
#$ -tc 5
#$ -N raccess_exact
#$ -o logs/uge/raccess_exact_$TASK_ID.out
#$ -e logs/uge/raccess_exact_$TASK_ID.err

set -euo pipefail

PROJECT_ROOT="/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper"
cd "${PROJECT_ROOT}"

BATCH_NUM=$(printf "%02d" ${SGE_TASK_ID})
BATCH_FILE="${PROJECT_ROOT}/benchmark/data/raccess_exact_batches/batch_${BATCH_NUM}.txt"
RESULTS_DIR="${PROJECT_ROOT}/benchmark/results/raccess_exact"
mkdir -p "${RESULTS_DIR}"
mkdir -p "logs/uge"

echo "Batch ${BATCH_NUM}: Processing sequences from ${BATCH_FILE}"

N_SEQS=$(wc -l < "${BATCH_FILE}")
echo "Total sequences in batch: ${N_SEQS}"

while IFS=$'\t' read -r SEQ_ID SEQ_LENGTH; do
    echo "Processing: ${SEQ_ID} (length=${SEQ_LENGTH})"
    
    TEMP_DIR=$(mktemp -d)
    ${PROJECT_ROOT}/.venv/bin/python3 << EOF
from Bio import SeqIO
for record in SeqIO.parse("${PROJECT_ROOT}/benchmark/data/bpRNA_1m_90.fasta", "fasta"):
    if record.id == "${SEQ_ID}":
        SeqIO.write([record], "${TEMP_DIR}/seq.fa", "fasta")
        break
EOF
    
    # Run Raccess-exact (max_span = sequence length)
    /usr/bin/time -v ${PROJECT_ROOT}/vendor/raccess/src/raccess/run_raccess \
      -seqfile="${TEMP_DIR}/seq.fa" \
      -access_len=7 \
      -max_span="${SEQ_LENGTH}" \
      -outfile="${TEMP_DIR}/out.txt" \
      2>"${TEMP_DIR}/time.log" || true
    
    # Parse results
    ${PROJECT_ROOT}/.venv/bin/python3 << PYEOF
import re
import json

result = {
    "seq_id": "${SEQ_ID}",
    "seq_length": ${SEQ_LENGTH},
    "method": "raccess-exact",
    "max_span": ${SEQ_LENGTH},
    "elapsed_sec": None,
    "memory_kb": None,
    "status": "FAILED",
    "exit_code": 1
}

with open("${TEMP_DIR}/time.log") as f:
    stderr = f.read()

if "Exit status: 0" in stderr:
    result["status"] = "SUCCESS"
    result["exit_code"] = 0

# Parse time (exact may take hours:minutes:seconds)
m = re.search(r'Elapsed.*?(\d+):(\d+):([\d.]+)', stderr)
if m:
    hours, minutes, seconds = int(m.group(1)), int(m.group(2)), float(m.group(3))
    result["elapsed_sec"] = hours * 3600 + minutes * 60 + seconds
else:
    m = re.search(r'Elapsed.*?(\d+):([\d.]+)', stderr)
    if m:
        minutes, seconds = int(m.group(1)), float(m.group(2))
        result["elapsed_sec"] = minutes * 60 + seconds

m = re.search(r'Maximum resident set size.*?(\d+)', stderr)
if m:
    result["memory_kb"] = int(m.group(1))

with open("${RESULTS_DIR}/${SEQ_ID}.json", "w") as f:
    json.dump(result, f, indent=2)
PYEOF
    
    rm -rf "${TEMP_DIR}"
    
done < "${BATCH_FILE}"

echo "Batch ${BATCH_NUM} completed"
