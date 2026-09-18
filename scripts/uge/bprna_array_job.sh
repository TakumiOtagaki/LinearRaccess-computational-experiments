#!/usr/bin/env bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=01:00:00
#$ -N linearraccess_bprna_array
#$ -t 1-26063
#$ -o logs/bprna_$JOB_ID.$TASK_ID.out
#$ -e logs/bprna_$JOB_ID.$TASK_ID.err
#$ -tc 10

set -euo pipefail

# Configuration
REPO_ROOT="/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper"
BPRNA_FASTA="$REPO_ROOT/benchmark/data/bpRNA_1m_90_filtered_clean.fa"
OUTDIR="$REPO_ROOT/benchmark/results/bprna_validation"
LOGDIR="$OUTDIR/logs"
LINRACC_BIN="$REPO_ROOT/vendor/LinearRaccess/build-raccess/LinRacc"

# Settings (production defaults)
BEAM=100
C_MULTI=30
C_HAIRPIN=30
ACCESS_LENS="1,3,7,10"
ENGINE="lincapr"
ENERGY="turner2004"

mkdir -p "$OUTDIR/raw" "$LOGDIR"

# Check prerequisites
if [[ ! -f "$BPRNA_FASTA" ]]; then
    echo "ERROR: bpRNA FASTA not found: $BPRNA_FASTA" >&2
    exit 1
fi

if [[ ! -x "$LINRACC_BIN" ]]; then
    echo "ERROR: LinRacc binary not found or not executable: $LINRACC_BIN" >&2
    exit 1
fi

# Get sequence for this task
TOTAL=$(grep -c "^>" "$BPRNA_FASTA" || true)
if [[ "$SGE_TASK_ID" -gt "$TOTAL" ]]; then
    echo "Task $SGE_TASK_ID exceeds total sequences ($TOTAL). Exiting." >&2
    exit 0
fi

# Extract sequence
SEQ_RECORD=$(awk -v n="$SGE_TASK_ID" '/^>/{i++} i==n{print; if(!/^>/) next}' "$BPRNA_FASTA")
SEQ_ID=$(echo "$SEQ_RECORD" | head -1 | sed 's/^>//' | awk '{print $1}')
SEQUENCE=$(echo "$SEQ_RECORD" | tail -n +2 | tr -d '\n' | tr 'T' 'U' | tr 'a-z' 'A-Z')

if [[ -z "$SEQ_ID" ]] || [[ -z "$SEQUENCE" ]]; then
    echo "ERROR: Failed to extract sequence for task $SGE_TASK_ID" >&2
    exit 1
fi

SEQ_LEN=${#SEQUENCE}
SAFE_ID=$(echo "$SEQ_ID" | tr '/ ' '_')

# Output paths
RAW_OUT="$OUTDIR/raw/${SAFE_ID}.txt"
JOB_LOG="$LOGDIR/${SAFE_ID}.log"
TMP_INPUT="$OUTDIR/raw/.tmp_${SAFE_ID}_$$.fa"

# Cleanup function
cleanup() {
    rm -f "$TMP_INPUT"
}
trap cleanup EXIT

# Create temporary input file
echo ">$SEQ_ID" > "$TMP_INPUT"
echo "$SEQUENCE" >> "$TMP_INPUT"

# Run LinearRaccess
echo "=== LinearRaccess Task $SGE_TASK_ID ===" >&2
echo "Date: $(date)" >&2
echo "Sequence ID: $SEQ_ID" >&2
echo "Length: $SEQ_LEN nt" >&2
echo "Settings: beam=$BEAM, c_multi=$C_MULTI, c_hairpin=$C_HAIRPIN" >&2

start_time=$(date +%s)
set +e
"$LINRACC_BIN" \
    -seqfile="$TMP_INPUT" \
    -outfile="$RAW_OUT" \
    -access_len="$ACCESS_LENS" \
    -engine="$ENGINE" \
    -beam="$BEAM" \
    -c_multi="$C_MULTI" \
    -c_hairpin="$C_HAIRPIN" \
    -energy="$ENERGY" \
    -probabilities \
    -metadata \
    > "$JOB_LOG" 2>&1
STATUS=$?
set -e
end_time=$(date +%s)

ELAPSED=$((end_time - start_time))

# Record result
RESULT_FILE="$OUTDIR/result_summary.tsv"
if [[ "$SGE_TASK_ID" -eq 1 ]]; then
    echo -e "task_id\tseq_id\tlength\tstatus\telapsed_sec\toutput_file" > "$RESULT_FILE"
fi

if [[ $STATUS -eq 0 ]]; then
    echo -e "$SGE_TASK_ID\t$SEQ_ID\t$SEQ_LEN\tSUCCESS\t$ELAPSED\t$RAW_OUT" >> "$RESULT_FILE"
    echo "SUCCESS: Completed in $ELAPSED seconds" >&2
else
    echo -e "$SGE_TASK_ID\t$SEQ_ID\t$SEQ_LEN\tFAILED\t$ELAPSED\t$RAW_OUT" >> "$RESULT_FILE"
    echo "ERROR: LinRacc failed with exit code $STATUS" >&2
    cat "$JOB_LOG" >&2
    exit $STATUS
fi

exit 0
