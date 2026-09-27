#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=06:00:00
#$ -t 1-8:1
#$ -N bprna_acc_capb
#$ -o logs/uge/bprna_acc_capb_$JOB_ID.$TASK_ID.out
#$ -j y
#
# Table 1 / Fig. 2 / Fig. S1-S2 accuracy benchmark, rerun for LinearRaccess
# with loop caps tied to the beam (c_multi = c_hairpin = b, the new default)
# and with w=35 added (the window length of DeepRaccess's published
# checkpoint).  Same protocol as bprna_accuracy_comparison.sh -- 12,872
# bpRNA sequences of 100-1,000 nt, Turner 1999, Raccess Exact as reference --
# but all four methods run with w in {1,3,7,10,35}, into a new results
# directory:
#   exact/        Raccess Exact (max_span = sequence length)
#   window_w100/  Raccess Window (max_span = 100)
#   linear_b100/  LinearRaccess, -engine=raccess, beam=100, caps auto
#   linear_b200/  LinearRaccess, -engine=raccess, beam=200, caps auto
# RNAplfold (-u 35) runs separately: rnaplfold_comparison.sh with ULENGTH=35.
#
# The LinRacc binary is frozen in ${RESULTS_DIR}/bin before submission (see
# the bottom of this file); per-sequence *_summary.json files make the job
# resumable.  Per-sequence timing here is incidental (16 runs share a node);
# runtime figures come from the single-core runtime jobs, not from this one.

set -euo pipefail

PROJECT_ROOT=/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper
cd "${PROJECT_ROOT}"

FASTA="${PROJECT_ROOT}/benchmark/data/bprna_100_1000.fa"
SEQ_LIST="${PROJECT_ROOT}/benchmark/data/bprna_100_1000_list.txt"
RESULTS_DIR="${PROJECT_ROOT}/benchmark/results/bprna_accuracy_capb"
RACCESS_BIN="${PROJECT_ROOT}/vendor/raccess/src/raccess/run_raccess"
LINRACC_BIN="${RESULTS_DIR}/bin/LinRacc"
ACCESS_LENS="1,3,7,10,35"
BEAMS="100 200"
N_PARALLEL=16
NUM_TASKS=8

test -x "${LINRACC_BIN}" || { echo "missing frozen binary ${LINRACC_BIN}" >&2; exit 1; }
mkdir -p "${RESULTS_DIR}"/{exact,window_w100,linear_b100,linear_b200,.work}

TOTAL_SEQ=$(wc -l < "${SEQ_LIST}")
TASK_ID="${SGE_TASK_ID:-1}"
[ "${TASK_ID}" = "undefined" ] && TASK_ID=1
BATCH=$(( (TOTAL_SEQ + NUM_TASKS - 1) / NUM_TASKS ))
START_IDX=$(( (TASK_ID - 1) * BATCH + 1 ))
END_IDX=$(( TASK_ID * BATCH < TOTAL_SEQ ? TASK_ID * BATCH : TOTAL_SEQ ))
echo "Start: $(date -u +%FT%TZ) host=$(hostname) task=${TASK_ID} seqs ${START_IDX}-${END_IDX}"

# Runs "$@" under /usr/bin/time -v; prints "<exit> <max_rss_kb> <wall_s>".
timed() {
    local log=$1; shift
    local rc=0
    /usr/bin/time -v "$@" > /dev/null 2> "${log}" || rc=$?
    local rss wall
    rss=$(awk -F': ' '/Maximum resident set size/ {print $2}' "${log}")
    wall=$(awk -F': ' '/Elapsed \(wall clock\)/ {print $2}' "${log}" \
           | awk -F: '{ if (NF == 3) print $1*3600 + $2*60 + $3; else print $1*60 + $2 }')
    echo "${rc} ${rss:-null} ${wall:-null}"
}

process_sequence() {
    local line_num=$1
    local seq_id seq_len
    read -r seq_id seq_len < <(sed -n "${line_num}p" "${SEQ_LIST}")
    [ -f "${RESULTS_DIR}/${seq_id}_summary.json" ] && return 0

    local work="${RESULTS_DIR}/.work/${seq_id}"
    mkdir -p "${work}"
    local seq_file="${work}/seq.fa"
    awk -v id=">${seq_id}" '
        /^>/ { keep = ($1 == id); if (keep) print $1; next }
        keep { printf "%s", $0 }
        END { print "" }' "${FASTA}" > "${seq_file}"

    local exact window b100 b200
    exact=$(timed "${work}/exact.log" "${RACCESS_BIN}" -seqfile="${seq_file}" \
        -access_len="${ACCESS_LENS}" -max_span="${seq_len}" \
        -outfile="${RESULTS_DIR}/exact/${seq_id}.txt")
    window=$(timed "${work}/window.log" "${RACCESS_BIN}" -seqfile="${seq_file}" \
        -access_len="${ACCESS_LENS}" -max_span=100 \
        -outfile="${RESULTS_DIR}/window_w100/${seq_id}.txt")
    local linear_json=""
    for beam in ${BEAMS}; do
        local r
        r=$(timed "${work}/linear_b${beam}.log" "${LINRACC_BIN}" -seqfile="${seq_file}" \
            -outfile="${RESULTS_DIR}/linear_b${beam}/${seq_id}.txt" \
            -access_len="${ACCESS_LENS}" -engine=raccess -energy=turner1999 \
            -beam="${beam}" -c_multi=auto -c_hairpin=auto -probabilities -metadata)
        read -r rc rss wall <<< "${r}"
        linear_json+=$(printf ',"linearraccess-b%s":{"exit_code":%s,"memory_kb":%s,"wall_s":%s,"beam":%s,"caps":"auto"}' \
                       "${beam}" "${rc}" "${rss}" "${wall}" "${beam}")
    done

    read -r e_rc e_rss e_wall <<< "${exact}"
    read -r w_rc w_rss w_wall <<< "${window}"
    printf '{"seq_id":"%s","seq_length":%s,"access_lens":"%s","methods":{"raccess-exact":{"exit_code":%s,"memory_kb":%s,"wall_s":%s,"max_span":%s},"raccess-w100":{"exit_code":%s,"memory_kb":%s,"wall_s":%s,"max_span":100}%s}}\n' \
        "${seq_id}" "${seq_len}" "${ACCESS_LENS}" "${e_rc}" "${e_rss}" "${e_wall}" "${seq_len}" \
        "${w_rc}" "${w_rss}" "${w_wall}" "${linear_json}" \
        > "${RESULTS_DIR}/${seq_id}_summary.json.tmp"
    mv "${RESULTS_DIR}/${seq_id}_summary.json.tmp" "${RESULTS_DIR}/${seq_id}_summary.json"
    rm -rf "${work}"
}

export -f process_sequence timed
export SEQ_LIST FASTA RESULTS_DIR RACCESS_BIN LINRACC_BIN ACCESS_LENS BEAMS

seq "${START_IDX}" "${END_IDX}" | xargs -P "${N_PARALLEL}" -I{} bash -c 'process_sequence {}'

n_done=$(ls "${RESULTS_DIR}" | grep -c '_summary.json$' || true)
echo "End: $(date -u +%FT%TZ) task=${TASK_ID} summaries_total=${n_done}"

# Submission (from the repository root):
#   mkdir -p benchmark/results/bprna_accuracy_capb/bin
#   cp vendor/LinearRaccess/build-raccess/LinRacc benchmark/results/bprna_accuracy_capb/bin/
#   LinRacc --build-info > benchmark/results/bprna_accuracy_capb/bin/build_info.json
#   qsub benchmark/scripts/uge/bprna_accuracy_capb.sh
