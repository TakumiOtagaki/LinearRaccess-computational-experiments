#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=8G
#$ -l mem_req=8G
#$ -t 1-25300:1
#$ -tc 10
#$ -N raccess_exact
#$ -o logs/uge/raccess_exact_$TASK_ID.out
#$ -e logs/uge/raccess_exact_$TASK_ID.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${PROJECT_ROOT}"

RESULTS_DIR="benchmark/results/raccess_exact"
mkdir -p "${RESULTS_DIR}"
mkdir -p "logs/uge"

# Get sequences <= 1000 nt from manifest
SEQ_LIST=$(python3 -c "
import json
with open('benchmark/results/bprna_runs/8666393/manifest.json') as f:
    data = json.load(f)
seqs = [s for s in data['sequences'] if s['length'] <= 1000]
print('\n'.join([s['id'] for s in seqs]))
")

# Get sequence for this task
SEQ_ID=$(echo "${SEQ_LIST}" | sed -n "${SGE_TASK_ID}p")

if [ -z "${SEQ_ID}" ]; then
    echo "ERROR: No sequence for task ${SGE_TASK_ID}"
    exit 1
fi

SEQ_LENGTH=$(python3 -c "
import json
with open('benchmark/results/bprna_runs/8666393/manifest.json') as f:
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

# Run Raccess-exact (max_span = sequence length)
echo "Running Raccess-exact (max_span=${SEQ_LENGTH})..."
/usr/bin/time -v vendor/raccess/src/raccess/run_raccess \
  -seqfile="${TEMP_DIR}/seq.fa" \
  -access_len=7 \
  -max_span="${SEQ_LENGTH}" \
  -outfile="${TEMP_DIR}/out.txt" \
  2>"${TEMP_DIR}/time.log" || true

# Parse results
uv run python3 << EOF
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

# Parse time log
with open("${TEMP_DIR}/time.log") as f:
    stderr = f.read()

# Check exit status
if "Exit status: 0" in stderr:
    result["status"] = "SUCCESS"
    result["exit_code"] = 0

# Parse elapsed time (exact may take hours)
m = re.search(r'Elapsed.*?(\d+):(\d+):([\d.]+)', stderr)
if m:
    hours, minutes, seconds = int(m.group(1)), int(m.group(2)), float(m.group(3))
    result["elapsed_sec"] = hours * 3600 + minutes * 60 + seconds
else:
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
