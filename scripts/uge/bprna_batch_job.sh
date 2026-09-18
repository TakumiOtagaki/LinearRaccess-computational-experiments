#!/usr/bin/env bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=02:00:00
#$ -N linearraccess_bprna_batch
#$ -t 1-27
#$ -o logs/batch_$JOB_ID.$TASK_ID.out
#$ -e logs/batch_$JOB_ID.$TASK_ID.err
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

# Batch processing: each job handles 1000 sequences
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

echo "=== Batch $SGE_TASK_ID: Processing sequences $START_IDX to $END_IDX ===" >&2

# Process each sequence in the batch
SUCCESS_COUNT=0
FAILED_COUNT=0

for IDX in $(seq $START_IDX $END_IDX); do
    # Extract sequence
    SEQ_RECORD=$(awk -v n="$IDX" '/^>/{i++} i==n{print; if(!/^>/) next}' "$BPRNA_FASTA")
    SEQ_ID=$(echo "$SEQ_RECORD" | head -1 | sed 's/^>//' | awk '{print $1}')
    SEQUENCE=$(echo "$SEQ_RECORD" | tail -n +2 | tr -d '\n')
    
    if [[ -z "$SEQ_ID" ]] || [[ -z "$SEQUENCE" ]]; then
        echo "WARNING: Failed to extract sequence $IDX" >&2
        ((FAILED_COUNT++))
        continue
    fi
    
    SEQ_LEN=${#SEQUENCE}
    SAFE_ID=$(echo "$SEQ_ID" | tr '/ ' '_')
    
    # Output paths
    RAW_OUT="$OUTDIR/raw/${SAFE_ID}.txt"
    JOB_LOG="$LOGDIR/${SAFE_ID}.log"
    TMP_INPUT="$OUTDIR/raw/.tmp_${SAFE_ID}_$$.fa"
    
    # Create temporary input file
    echo ">$SEQ_ID" > "$TMP_INPUT"
    echo "$SEQUENCE" >> "$TMP_INPUT"
    
    # Run LinearRaccess
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
    
    # Cleanup
    rm -f "$TMP_INPUT"
    
    # Record result
    RESULT_FILE="$OUTDIR/batch_${SGE_TASK_ID}_summary.tsv"
    if [[ $IDX -eq $START_IDX ]]; then
        echo -e "batch_id\tseq_idx\tseq_id\tlength\tstatus\toutput_file" > "$RESULT_FILE"
    fi
    
    if [[ $STATUS -eq 0 ]]; then
        echo -e "$SGE_TASK_ID\t$IDX\t$SEQ_ID\t$SEQ_LEN\tSUCCESS\t$RAW_OUT" >> "$RESULT_FILE"
        ((SUCCESS_COUNT++))
    else
        echo -e "$SGE_TASK_ID\t$IDX\t$SEQ_ID\t$SEQ_LEN\tFAILED\t$RAW_OUT" >> "$RESULT_FILE"
        ((FAILED_COUNT++))
        echo "ERROR: Failed for $SEQ_ID (exit $STATUS)" >&2
    fi
done

echo "=== Batch $SGE_TASK_ID Complete ===" >&2
echo "Success: $SUCCESS_COUNT, Failed: $FAILED_COUNT" >&2

exit 0
