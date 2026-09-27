#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=08:00:00
#$ -t 1-3:1
#$ -N bprna_beam_capb
#$ -o logs/uge/bprna_beam_capb_$JOB_ID.$TASK_ID.out
#$ -j y
#
# Supplementary S3 beam-width sweep, rerun with loop caps tied to the beam
# (c_multi = c_hairpin = b) and w=35 added.  Same protocol as
# bprna_beam_sweep.sh: LinearRaccess only (-engine=raccess, Turner 1999) on
# the 1,770 bpRNA sequences of 300-1,000 nt, scored against Raccess Exact.
# Task t runs beam BEAMS[t]; beam=100 and beam=200 come from
# bprna_accuracy_capb.sh, and the exact reference from its exact/ directory.
# Uses the same frozen binary as bprna_accuracy_capb.sh.

set -euo pipefail

PROJECT_ROOT=/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper
cd "${PROJECT_ROOT}"

BEAMS=(50 300 500)
TASK_ID="${SGE_TASK_ID:-1}"
[ "${TASK_ID}" = "undefined" ] && TASK_ID=1
BEAM="${BEAMS[$((TASK_ID - 1))]}"

FASTA="${PROJECT_ROOT}/benchmark/data/bprna_100_1000.fa"
SEQ_LIST="${PROJECT_ROOT}/benchmark/data/bprna_300_1000_list.txt"
LINRACC_BIN="${PROJECT_ROOT}/benchmark/results/bprna_accuracy_capb/bin/LinRacc"
RESULTS_DIR="${PROJECT_ROOT}/benchmark/results/bprna_beam_sweep_300_1000_capb/linear_b${BEAM}"
ACCESS_LENS="1,3,7,10,35"
N_PARALLEL=16

test -x "${LINRACC_BIN}" || { echo "missing frozen binary ${LINRACC_BIN}" >&2; exit 1; }
mkdir -p "${RESULTS_DIR}/.work"
echo "Start: $(date -u +%FT%TZ) host=$(hostname) beam=${BEAM}"

process_sequence() {
    local line_num=$1
    local seq_id seq_len
    read -r seq_id seq_len < <(sed -n "${line_num}p" "${SEQ_LIST}")
    [ -f "${RESULTS_DIR}/${seq_id}.txt" ] && return 0
    local seq_file="${RESULTS_DIR}/.work/${seq_id}.fa"
    awk -v id=">${seq_id}" '
        /^>/ { keep = ($1 == id); if (keep) print $1; next }
        keep { printf "%s", $0 }
        END { print "" }' "${FASTA}" > "${seq_file}"
    "${LINRACC_BIN}" -seqfile="${seq_file}" -outfile="${RESULTS_DIR}/.work/${seq_id}.txt" \
        -access_len="${ACCESS_LENS}" -engine=raccess -energy=turner1999 \
        -beam="${BEAM}" -c_multi=auto -c_hairpin=auto -probabilities -metadata \
        > /dev/null 2>&1
    mv "${RESULTS_DIR}/.work/${seq_id}.txt" "${RESULTS_DIR}/${seq_id}.txt"
    rm -f "${seq_file}"
}

export -f process_sequence
export SEQ_LIST FASTA RESULTS_DIR LINRACC_BIN ACCESS_LENS BEAM

seq 1 "$(wc -l < "${SEQ_LIST}")" | xargs -P "${N_PARALLEL}" -I{} bash -c 'process_sequence {}'

echo "End: $(date -u +%FT%TZ) beam=${BEAM} outputs=$(ls "${RESULTS_DIR}" | grep -c '\.txt$')"
