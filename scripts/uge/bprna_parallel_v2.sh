#!/usr/bin/env bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=04:00:00
#$ -N linearraccess_bprna_v2
#$ -t 1-27
#$ -o logs/v2_$JOB_ID.$TASK_ID.out
#$ -e logs/v2_$JOB_ID.$TASK_ID.err
#$ -tc 10

set -euo pipefail

# Configuration
REPO_ROOT="/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper"
BPRNA_FASTA="$REPO_ROOT/benchmark/data/bpRNA_1m_90_filtered_clean.fa"
OUTDIR="$REPO_ROOT/benchmark/results/bprna_parallel_v2"
LOGDIR="$OUTDIR/logs"
LINRACC_BIN="$REPO_ROOT/vendor/LinearRaccess/build-raccess/LinRacc"

# Settings (production defaults)
BEAM=100
C_MULTI=30
C_HAIRPIN=30
ACCESS_LENS="1,3,7,10"
ENGINE="lincapr"
ENERGY="turner2004"
PARALLEL=10

mkdir -p "$OUTDIR/raw" "$LOGDIR"

echo "=== Job $SGE_TASK_ID Starting ===" >&2

# Use seqkit to extract and process sequences in parallel if available
# Otherwise use awk with xargs
if command -v seqkit &> /dev/null; then
    echo "Using seqkit for sequence extraction" >&2
    seqkit split -p 10000 -O "$OUTDIR/tmp_split" "$BPRNA_FASTA" 2>/dev/null || true
fi

# Alternative: Pre-split the FASTA into individual files for this batch
BATCH_START=$(( (SGE_TASK_ID - 1) * 1000 + 1 ))
BATCH_END=$(( SGE_TASK_ID * 1000 ))
TOTAL_SEQ=26063
if [[ $BATCH_END -gt $TOTAL_SEQ ]]; then
    BATCH_END=$TOTAL_SEQ
fi

echo "Processing sequences $BATCH_START to $BATCH_END" >&2

# Create output directory for this batch
BATCH_OUTDIR="$OUTDIR/batch_$SGE_TASK_ID"
mkdir -p "$BATCH_OUTDIR" "$BATCH_OUTDIR/logs"

# Initialize summary
echo -e "seq_idx\tseq_id\tlength\tstatus\telapsed_sec" > "$BATCH_OUTDIR/summary.tsv"

# Function to process one sequence
process_one() {
    local IDX=$1
    local OUTDIR=$2
    local LOGDIR=$3

    # Extract the IDX-th sequence using awk
    awk -v n="$IDX" '/^>/{i++} i==n{print; if(!/^>/) next}' "$BPRNA_FASTA" 2>/dev/null | \
    awk 'BEGIN{RS=">"; ORS=""} NR==2{print ">" $0}' > "$OUTDIR/tmp_$$_${IDX}.fa" 2>/dev/null

    if [[ ! -s "$OUTDIR/tmp_$$_${IDX}.fa" ]]; then
        rm -f "$OUTDIR/tmp_$$_${IDX}.fa"
        return 1
    fi

    local SEQ_ID=$(head -1 "$OUTDIR/tmp_$$_${IDX}.fa" | sed 's/^>//' | awk '{print $1}')
    local SEQ_LEN=$(tail -n +2 "$OUTDIR/tmp_$$_${IDX}.fa" | tr -d '\n' | wc -c)

    local START_TIME=$(date +%s)

    if "$LINRACC_BIN" \
        -seqfile="$OUTDIR/tmp_$$_${IDX}.fa" \
        -outfile="$OUTDIR/${SEQ_ID}.txt" \
        -access_len="$ACCESS_LENS" \
        -engine="$ENGINE" \
        -beam="$BEAM" \
        -c_multi="$C_MULTI" \
        -c_hairpin="$C_HAIRPIN" \
        -energy="$ENERGY" \
        -probabilities \
        -metadata \
        > "$LOGDIR/${SEQ_ID}.log" 2>&1; then
        local END_TIME=$(date +%s)
        local ELAPSED=$((END_TIME - START_TIME))
        echo -e "$IDX\t$SEQ_ID\t$SEQ_LEN\tSUCCESS\t$ELAPSED" >> "$BATCH_OUTDIR/summary.tsv"
        rm -f "$OUTDIR/tmp_$$_${IDX}.fa"
        return 0
    else
        local END_TIME=$(date +%s)
        local ELAPSED=$((END_TIME - START_TIME))
        echo -e "$IDX\t$SEQ_ID\t$SEQ_LEN\tFAILED\t$ELAPSED" >> "$BATCH_OUTDIR/summary.tsv"
        rm -f "$OUTDIR/tmp_$$_${IDX}.fa"
        return 1
    fi
}

export -f process_one
export BPRNA_FASTA LINRACC_BIN BEAM C_MULTI C_HAIRPIN ACCESS_LENS ENGINE ENERGY

# Use GNU parallel if available, otherwise use simple loop with background jobs
if command -v parallel &> /dev/null; then
    seq $BATCH_START $BATCH_END | parallel -j $PARALLEL process_one {} "$BATCH_OUTDIR" "$LOGDIR"
else
    # Simple background job approach
    ACTIVE_JOBS=0
    for IDX in $(seq $BATCH_START $BATCH_END); do
        while [[ $ACTIVE_JOBS -ge $PARALLEL ]]; do
            wait -n 2>/dev/null || true
            ACTIVE_JOBS=$((ACTIVE_JOBS - 1))
        done

        process_one "$IDX" "$BATCH_OUTDIR" "$LOGDIR" &
        ACTIVE_JOBS=$((ACTIVE_JOBS + 1))
    done
    wait
fi

# Aggregate results
SUCCESS=$(grep -c SUCCESS "$BATCH_OUTDIR/summary.tsv" 2>/dev/null || echo 0)
FAILED=$(grep -c FAILED "$BATCH_OUTDIR/summary.tsv" 2>/dev/null || echo 0)

echo "=== Job $SGE_TASK_ID Complete: $SUCCESS success, $FAILED failed ===" >&2

exit 0
