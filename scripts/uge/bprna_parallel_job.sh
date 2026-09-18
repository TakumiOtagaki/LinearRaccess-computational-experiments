#!/usr/bin/env bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=04:00:00
#$ -N linearraccess_bprna_parallel
#$ -t 1-27
#$ -o logs/parallel_$JOB_ID.$TASK_ID.out
#$ -e logs/parallel_$JOB_ID.$TASK_ID.err
#$ -tc 10

set -euo pipefail

# Configuration
REPO_ROOT="/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper"
BPRNA_FASTA="$REPO_ROOT/benchmark/data/bpRNA_1m_90_filtered_clean.fa"
OUTDIR="$REPO_ROOT/benchmark/results/bprna_validation_parallel"
LOGDIR="$OUTDIR/logs"
LINRACC_BIN="$REPO_ROOT/vendor/LinearRaccess/build-raccess/LinRacc"

# Settings (production defaults)
BEAM=100
C_MULTI=30
C_HAIRPIN=30
ACCESS_LENS="1,3,7,10"
ENGINE="lincapr"
ENERGY="turner2004"

# Parallel processing settings
PARALLEL_JOBS=10
BATCH_SIZE=1000
START_IDX=$(( (SGE_TASK_ID - 1) * BATCH_SIZE + 1 ))
END_IDX=$(( SGE_TASK_ID * BATCH_SIZE ))

# For the last batch, handle remaining sequences
TOTAL_SEQ=26063
if [[ $END_IDX -gt $TOTAL_SEQ ]]; then
    END_IDX=$TOTAL_SEQ
fi

mkdir -p "$OUTDIR/raw" "$LOGDIR"

# Check prerequisites
if [[ ! -f "$BPRNA_FASTA" ]]; then
    echo "ERROR: bpRNA FASTA not found: $BPRNA_FASTA" >&2
    exit 1
fi

if [[ ! -x "$LINRACC_BIN" ]]; then
    echo "ERROR: LinRacc binary not found: $LINRACC_BIN" >&2
    exit 1
fi

echo "=== Batch $SGE_TASK_ID: Processing sequences $START_IDX to $END_IDX (10 parallel) ===" >&2

# Create a temporary job list for this batch
JOB_LIST="$OUTDIR/.joblist_batch_${SGE_TASK_ID}_$$"
trap "rm -f '$JOB_LIST'" EXIT

# Extract sequences for this batch and create job list
awk -v start="$START_IDX" -v end="$END_IDX" '
    /^>/{n++; if(n>=start && n<=end){header=substr($0,2); id=$1; gsub(/^>/, "", id); print id}}
' "$BPRNA_FASTA" > "$JOB_LIST"

TOTAL_IN_BATCH=$((END_IDX - START_IDX + 1))
echo "Total sequences in batch: $TOTAL_IN_BATCH" >&2

# Function to process a single sequence
process_sequence() {
    local IDX=$1

    # Extract sequence
    local SEQ_RECORD=$(awk -v n="$IDX" '/^>/{i++} i==n{print; for(j=1;j<=10000;j++){getline; if(/^>/ || NF==0) exit; print}}' "$BPRNA_FASTA" 2>/dev/null)
    local SEQ_ID=$(echo "$SEQ_RECORD" | head -1 | sed 's/^>//' | awk '{print $1}')
    local SEQUENCE=$(echo "$SEQ_RECORD" | tail -n +2 | tr -d '\n' | tr -d ' ')

    if [[ -z "$SEQ_ID" ]] || [[ -z "$SEQUENCE" ]]; then
        echo "ERROR: Failed to extract sequence $IDX" >&2
        return 1
    fi

    local SEQ_LEN=${#SEQUENCE}
    local SAFE_ID=$(echo "$SEQ_ID" | tr '/ ' '_')

    # Output paths
    local RAW_OUT="$OUTDIR/raw/${SAFE_ID}.txt"
    local JOB_LOG="$LOGDIR/${SAFE_ID}.log"
    local TMP_INPUT="$OUTDIR/raw/.tmp_${SAFE_ID}_$IDX.fa"

    # Create temporary input file
    echo ">$SEQ_ID" > "$TMP_INPUT"
    echo "$SEQUENCE" >> "$TMP_INPUT"

    # Run LinearRaccess
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
    local STATUS=$?

    # Cleanup
    rm -f "$TMP_INPUT"

    # Record result
    if [[ $STATUS -eq 0 ]]; then
        echo -e "$SGE_TASK_ID\t$IDX\t$SEQ_ID\t$SEQ_LEN\tSUCCESS\t$RAW_OUT" >> "$OUTDIR/batch_${SGE_TASK_ID}_summary.tsv"
        echo "SUCCESS: $SEQ_ID ($SEQ_LEN nt)" >&2
    else
        echo -e "$SGE_TASK_ID\t$IDX\t$SEQ_ID\t$SEQ_LEN\tFAILED\t$RAW_OUT" >> "$OUTDIR/batch_${SGE_TASK_ID}_summary.tsv"
        echo "FAILED: $SEQ_ID (exit $STATUS)" >&2
    fi

    return $STATUS
}

export -f process_sequence
export BPRNA_FASTA OUTDIR LOGDIR LINRACC_BIN BEAM C_MULTI C_HAIRPIN ACCESS_LENS ENGINE ENERGY

# Initialize summary file
echo -e "batch_id\tseq_idx\tseq_id\tlength\tstatus\toutput_file" > "$OUTDIR/batch_${SGE_TASK_ID}_summary.tsv"

# Process sequences in parallel using background jobs
SUCCESS_COUNT=0
FAILED_COUNT=0
RUNNING_JOBS=0

for IDX in $(seq $START_IDX $END_IDX); do
    # Wait if we have too many parallel jobs
    while [[ $RUNNING_JOBS -ge $PARALLEL_JOBS ]]; do
        wait -n 2>/dev/null || true
        RUNNING_JOBS=$((RUNNING_JOBS - 1))
    done

    # Start background job
    (
        process_sequence "$IDX"
        if [[ $? -eq 0 ]]; then
            echo "SUCCESS"
        else
            echo "FAILED"
        fi
    ) &

    RUNNING_JOBS=$((RUNNING_JOBS + 1))
done

# Wait for all background jobs to complete
wait

# Count results
SUCCESS_COUNT=$(grep -c "SUCCESS" "$OUTDIR/batch_${SGE_TASK_ID}_summary.tsv" 2>/dev/null || echo 0)
FAILED_COUNT=$(grep -c "FAILED" "$OUTDIR/batch_${SGE_TASK_ID}_summary.tsv" 2>/dev/null || echo 0)

echo "=== Batch $SGE_TASK_ID Complete ===" >&2
echo "Success: $SUCCESS_COUNT, Failed: $FAILED_COUNT" >&2

exit 0
