#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=08:00:00
#$ -N linracc_bprna_acc
#$ -o logs/uge/linracc_bprna_accuracy_$JOB_ID.out
#$ -j y
#
# LinearRaccess-only accuracy run on a bpRNA list (Raccess outputs reused):
# -engine=raccess, Turner 1999, caps auto, w in {1,3,7,10,35}, beam $BEAM,
# into $OUTDIR/linear_b$BEAM, with the frozen binary $OUTDIR/bin/LinRacc.
# Required env: BEAM, LIST, OUTDIR (qsub -v).  16 concurrent runs, resumable.
set -euo pipefail
PROJECT_ROOT=/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper
cd "${PROJECT_ROOT}"
: "${BEAM:?}" "${LIST:?}" "${OUTDIR:?}"
FASTA=benchmark/data/bprna_100_1000.fa
BIN="${OUTDIR}/bin/LinRacc"
RES="${OUTDIR}/linear_b${BEAM}"
test -x "${BIN}"
mkdir -p "${RES}/.work"
echo "Start: $(date -u +%FT%TZ) host=$(hostname) beam=${BEAM} list=${LIST} $("${BIN}" --build-info)"
run_one() {
    local id=$1
    [ -f "${RES}/${id}.txt" ] && return 0
    awk -v id=">${id}" '/^>/ {keep = ($1 == id); if (keep) print $1; next} keep {printf "%s", $0} END {print ""}' "${FASTA}" > "${RES}/.work/${id}.fa"
    "${BIN}" -seqfile="${RES}/.work/${id}.fa" -outfile="${RES}/.work/${id}.txt" -access_len=1,3,7,10,35 \
        -engine=raccess -energy=turner1999 -beam="${BEAM}" -c_multi=auto -c_hairpin=auto \
        -probabilities -metadata > /dev/null 2>&1
    mv "${RES}/.work/${id}.txt" "${RES}/${id}.txt"
    rm -f "${RES}/.work/${id}.fa"
}
export -f run_one
export RES BIN BEAM FASTA
cut -f1 "${LIST}" | xargs -P 16 -I{} bash -c 'run_one {}'
echo "End: $(date -u +%FT%TZ) outputs=$(ls "${RES}" | grep -c '\.txt$')/$(wc -l < "${LIST}")"
