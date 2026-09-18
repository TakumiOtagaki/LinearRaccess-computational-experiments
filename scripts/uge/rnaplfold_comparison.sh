#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_4=1
#$ -l h_rt=02:00:00
#$ -t 1-16:1
#$ -tc 16
#$ -N rnaplfold_cmp
#$ -o logs/uge/rnaplfold_$TASK_ID.out
#$ -e logs/uge/rnaplfold_$TASK_ID.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${PROJECT_ROOT}"

# Same TMPDIR avoidance as bprna_accuracy_comparison.sh: every scratch path
# here is an explicit directory under the repo, nothing depends on TMPDIR.

if [ -z "${RNAPLFOLD_BIN:-}" ]; then
    module load ViennaRNA/2.7.2
    RNAPLFOLD_BIN="$(command -v RNAplfold)"
fi
PARAMETER_FILE="${RNAPLFOLD_PARAMETER_FILE:-$(dirname "${RNAPLFOLD_BIN}")/../share/ViennaRNA/rna_turner1999.par}"
test -r "${PARAMETER_FILE}"

FASTA="${PROJECT_ROOT}/benchmark/data/bprna_100_1000.fa"
RESULTS_DIR="${PROJECT_ROOT}/benchmark/results/rnaplfold_turner1999_w100_u10"
WORK_ROOT="${RESULTS_DIR}/.work"
LOGDIR="${PROJECT_ROOT}/logs/uge"

WINSIZE=100
SPAN=100
ULENGTH=10
TOTAL_SEQ=12872
NUM_TASKS=16

mkdir -p "${RESULTS_DIR}/lunp" "${RESULTS_DIR}/tasks" "${WORK_ROOT}" "${LOGDIR}"

if [ -n "${SGE_TASK_ID:-}" ] && [ "${SGE_TASK_ID}" != "undefined" ]; then
    TASK_ID="${SGE_TASK_ID}"
else
    TASK_ID="${TASK_ID:-1}"
fi
if ! [[ "${TASK_ID}" =~ ^[0-9]+$ ]] || (( TASK_ID < 1 || TASK_ID > NUM_TASKS )); then
    echo "Invalid task ID: ${TASK_ID}" >&2
    exit 1
fi
if [ -e "${RESULTS_DIR}/tasks/${TASK_ID}.json" ]; then
    echo "Task ${TASK_ID} already completed; refusing to overwrite frozen output" >&2
    exit 1
fi

BATCH_SIZE=$(( (TOTAL_SEQ + NUM_TASKS - 1) / NUM_TASKS ))
START_REC=$(( (TASK_ID - 1) * BATCH_SIZE + 1 ))
END_REC=$(( TASK_ID * BATCH_SIZE ))
if [ ${END_REC} -gt ${TOTAL_SEQ} ]; then
    END_REC=${TOTAL_SEQ}
fi

TASK_DIR="${WORK_ROOT}/task_${TASK_ID}"
mkdir -p "${TASK_DIR}"
CHUNK_FASTA="${TASK_DIR}/chunk.fa"

echo "Task ${TASK_ID}: records ${START_REC}-${END_REC} of ${TOTAL_SEQ}"

# Extract this task's slice of FASTA records (1-indexed, inclusive) using awk
# on a two-line-per-record file (the subset FASTA is single-line sequences).
awk -v start="${START_REC}" -v end="${END_REC}" '
    BEGIN { rec = 0 }
    /^>/ { rec++ }
    { if (rec >= start && rec <= end) print }
' "${FASTA}" > "${CHUNK_FASTA}"

N_IN_CHUNK=$(grep -c "^>" "${CHUNK_FASTA}")
echo "Extracted ${N_IN_CHUNK} sequences into ${CHUNK_FASTA}"

cd "${TASK_DIR}"
"${RNAPLFOLD_BIN}" -P "${PARAMETER_FILE}" -W "${WINSIZE}" -L "${SPAN}" -u "${ULENGTH}" < "${CHUNK_FASTA}"

# Move only the _lunp files (accessibility profiles) to the permanent
# results directory; the _dp.ps postscript dot plots are not needed and are
# left behind / cleaned up with the task directory.
python3 - "${CHUNK_FASTA}" "${RESULTS_DIR}" "${TASK_ID}" "${RNAPLFOLD_BIN}" "${PARAMETER_FILE}" <<'PY'
import hashlib
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

chunk, results, task, binary, parameters = sys.argv[1:]
chunk, results = Path(chunk), Path(results)
ids = [line[1:].split()[0] for line in chunk.read_text().splitlines() if line.startswith('>')]
sources = [chunk.parent / f'{seq_id}_lunp' for seq_id in ids]
if len(ids) != len(set(ids)) or any(not p.is_file() or not p.stat().st_size for p in sources):
    raise SystemExit('Duplicate sequence IDs or missing/empty RNAplfold output')
if any((results / 'lunp' / p.name).exists() for p in sources):
    raise SystemExit('Refusing to overwrite existing RNAplfold output')
for source in sources:
    shutil.move(str(source), results / 'lunp' / source.name)
record = {
    'task_id': int(task), 'n_sequences': len(ids),
    'completed_at': datetime.now(timezone.utc).isoformat(),
    'chunk_sha256': hashlib.sha256(chunk.read_bytes()).hexdigest(),
    'binary_sha256': hashlib.sha256(Path(binary).read_bytes()).hexdigest(),
    'parameter_sha256': hashlib.sha256(Path(parameters).read_bytes()).hexdigest(),
    'command': [binary, '-P', str(Path(parameters).resolve()), '-W', '100', '-L', '100', '-u', '10'],
}
(results / 'tasks' / f'{task}.json').write_text(json.dumps(record, indent=2) + '\n')
PY

echo "Task ${TASK_ID} done: ${N_IN_CHUNK} sequences processed"

cd "${PROJECT_ROOT}"
rm -rf "${TASK_DIR}"

exit 0
