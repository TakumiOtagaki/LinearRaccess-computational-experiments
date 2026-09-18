#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=8G
#$ -l mem_req=8G
#$ -t 1-160:1
#$ -tc 20
#$ -N raccess_cmp
#$ -o logs/uge/raccess_comparison_$TASK_ID.out
#$ -e logs/uge/raccess_comparison_$TASK_ID.err

set -euo pipefail

# Setup
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${PROJECT_ROOT}"

# Create output directory
RESULTS_DIR="benchmark/results/raccess_comparison"
mkdir -p "${RESULTS_DIR}"
mkdir -p "logs/uge"

# Get sequence ID for this task
SEQ_LIST="benchmark/data/raccess_comparison_sample.txt"
SEQ_ID=$(sed -n "${SGE_TASK_ID}p" "${SEQ_LIST}")

if [ -z "${SEQ_ID}" ]; then
    echo "ERROR: No sequence for task ${SGE_TASK_ID}"
    exit 1
fi

echo "Task ${SGE_TASK_ID}: Processing ${SEQ_ID}"

# Extract single sequence to temp file
TEMP_DIR=$(mktemp -d)
SEQ_FILE="${TEMP_DIR}/seq.fa"
uv run python3 << EOF
from Bio import SeqIO
for record in SeqIO.parse("benchmark/data/raccess_comparison_sample.fa", "fasta"):
    if record.id == "${SEQ_ID}":
        SeqIO.write([record], "${SEQ_FILE}", "fasta")
        print(f"Length: {len(record.seq)}")
        break
EOF

SEQ_LENGTH=$(grep -v "^>" "${SEQ_FILE}" | tr -d '\n' | wc -c)
echo "Sequence length: ${SEQ_LENGTH}"

# Run Raccess with different max_span values
# 1. Raccess-exact (max_span = sequence length)
# 2. Raccess-W100
# 3. Raccess-W200

METHODS=("exact" "W100" "W200")
MAX_SPANS=(${SEQ_LENGTH} 100 200)

for i in ${!METHODS[@]}; do
    METHOD=${METHODS[$i]}
    MAX_SPAN=${MAX_SPANS[$i]}
    
    echo "Running Raccess-${METHOD} (max_span=${MAX_SPAN})..."
    
    # Skip exact for very long sequences (>1000 nt) as it will be too slow
    if [ "${METHOD}" = "exact" ] && [ ${SEQ_LENGTH} -gt 1000 ]; then
        echo "Skipping exact for long sequence (>1000 nt)"
        cat > "${RESULTS_DIR}/${SEQ_ID}_${METHOD}.json" << EOF
{
  "seq_id": "${SEQ_ID}",
  "seq_length": ${SEQ_LENGTH},
  "method": "raccess-${METHOD}",
  "access_len": "7",
  "max_span": ${MAX_SPAN},
  "elapsed_sec": null,
  "memory_kb": null,
  "status": "SKIPPED (too long)",
  "exit_code": -1
}
EOF
        continue
    fi
    
    uv run python3 benchmark/scripts/run_raccess_single.py \
        --seq-file "${SEQ_FILE}" \
        --seq-id "${SEQ_ID}" \
        --seq-length ${SEQ_LENGTH} \
        --access-len 7 \
        --max-span ${MAX_SPAN} \
        --output-dir "${RESULTS_DIR}" \
        --method "raccess-${METHOD}" \
        --output-json "${RESULTS_DIR}/${SEQ_ID}_${METHOD}.json" \
        || echo "Failed: ${SEQ_ID} ${METHOD}"
    
    echo "Completed Raccess-${METHOD}"
done

# Run LinearRaccess for comparison
echo "Running LinearRaccess..."
uv run python3 benchmark/scripts/run_linearraccess_single.py \
    --seq-file "${SEQ_FILE}" \
    --seq-id "${SEQ_ID}" \
    --seq-length ${SEQ_LENGTH} \
    --access-len 7 \
    --beam 200 \
    --output-dir "${RESULTS_DIR}" \
    --method "linearraccess-beam200" \
    --output-json "${RESULTS_DIR}/${SEQ_ID}_linearraccess.json" \
    || echo "Failed: ${SEQ_ID} linearraccess"

# Cleanup
rm -rf "${TEMP_DIR}"

echo "Task ${SGE_TASK_ID} completed for ${SEQ_ID}"
