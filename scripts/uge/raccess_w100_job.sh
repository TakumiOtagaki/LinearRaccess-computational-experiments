#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=4G
#$ -l mem_req=4G
#$ -t 1-26063:1
#$ -tc 20
#$ -N raccess_w100
#$ -o logs/uge/raccess_w100_$TASK_ID.out
#$ -e logs/uge/raccess_w100_$TASK_ID.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${PROJECT_ROOT}"

RESULTS_DIR="benchmark/results/raccess_w100"
mkdir -p "${RESULTS_DIR}"
mkdir -p "logs/uge"

# Get sequence from manifest
MANIFEST="benchmark/results/bprna_runs/8666393/manifest.json"
SEQ_ID=$(python3 -c "
import json
with open('${MANIFEST}') as f:
    data = json.load(f)
seqs = data['sequences']
idx = ${SGE_TASK_ID} - 1
if idx < len(seqs):
    print(seqs[idx]['id'])
")

if [ -z "${SEQ_ID}" ]; then
    echo "ERROR: No sequence for task ${SGE_TASK_ID}"
    exit 1
fi

SEQ_LENGTH=$(python3 -c "
import json
with open('${MANIFEST}') as f:
    data = json.load(f)
for s in data['sequences']:
    if s['id'] == '${SEQ_ID}':
        print(s['length'])
        break
")

echo "Task ${SGE_TASK_ID}: ${SEQ_ID} (length=${SEQ_LENGTH})"

# Extract sequence to temp file
TEMP_DIR=$(mktemp -d)
uv run python3 << EOF
from Bio import SeqIO
for record in SeqIO.parse("benchmark/data/bpRNA_1m_90.fasta", "fasta"):
    if record.id == "${SEQ_ID}":
        SeqIO.write([record], "${TEMP_DIR}/seq.fa", "fasta")
        break
EOF

# Run Raccess-W100
echo "Running Raccess-W100..."
/usr/bin/time -v vendor/raccess/src/raccess/run_raccess \
  -seqfile="${TEMP_DIR}/seq.fa" \
  -access_len=7 \
  -max_span=100 \
  -outfile="${TEMP_DIR}/out.txt" \
  2>"${TEMP_DIR}/time.log" || true

# Parse results
uv run python3 << EOF
import re
import json

result = {
    "seq_id": "${SEQ_ID}",
    "seq_length": ${SEQ_LENGTH},
    "method": "raccess-w100",
    "max_span": 100,
    "elapsed_sec": None,
    "memory_kb": None,
    "status": "FAILED",
    "exit_code": 1
}

# Parse time log
with open("${TEMP_DIR}/time.log") as f:
    stderr = f.read()

# Check exit status
if "Exit status: 0" in stderr:
    result["status"] = "SUCCESS"
    result["exit_code"] = 0

# Parse elapsed time
m = re.search(r'Elapsed.*?(\d+):([\d.]+)', stderr)
if m:
    minutes, seconds = int(m.group(1)), float(m.group(2))
    result["elapsed_sec"] = minutes * 60 + seconds

# Parse memory
m = re.search(r'Maximum resident set size.*?(\d+)', stderr)
if m:
    result["memory_kb"] = int(m.group(1))

# Save
with open("${RESULTS_DIR}/${SEQ_ID}.json", "w") as f:
    json.dump(result, f, indent=2)

print(json.dumps(result, indent=2))
EOF

# Cleanup
rm -rf "${TEMP_DIR}"

echo "Task ${SGE_TASK_ID} completed"
