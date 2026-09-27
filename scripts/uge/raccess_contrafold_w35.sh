#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=12:00:00
#$ -N racc_cf_w35
#$ -o logs/uge/raccess_contrafold_w35_$JOB_ID.out
#$ -j y
#
# DeepRaccess comparison references at w=35: legacy Raccess with the
# CONTRAfold model (run_raccess_contrafold), which is DeepRaccess's teacher.
#   w100/   -max_span=100 (the teacher setting for FCN_structured.pth; this
#           reproduces the shipped sample_data teacher CSVs exactly)
#   exact/  -max_span=N   (no span limit)
# for the 12,872 Table 1 bpRNA sequences and the 20 long (1,000-3,541 nt) ones.
set -euo pipefail
PROJECT_ROOT=/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper
cd "${PROJECT_ROOT}"
OUT=benchmark/results/deepraccess_comparison/raccess_contrafold
BIN=vendor/raccess/src/raccess/run_raccess_contrafold
mkdir -p "${OUT}"/{w100,exact,.work}
sha256sum "${BIN}" > "${OUT}/sha256.txt"
echo "Start: $(date -u +%FT%TZ) host=$(hostname)"

run_one() {  # <id> <len> <fasta>
    local id=$1 len=$2 fasta=$3
    [ -f "${OUT}/exact/${id}.txt" ] && return 0
    local fa="${OUT}/.work/${id}.fa"
    awk -v id=">${id}" '/^>/ {keep = ($1 == id); if (keep) print $1; next} keep {printf "%s", $0} END {print ""}' "${fasta}" > "${fa}"
    "${BIN}" -seqfile="${fa}" -access_len=35 -max_span=100 -outfile="${OUT}/w100/${id}.txt" > /dev/null
    /usr/bin/time -f "%e %M" -o "${OUT}/.work/${id}.time" \
        "${BIN}" -seqfile="${fa}" -access_len=35 -max_span="${len}" -outfile="${OUT}/.work/${id}.exact" > /dev/null
    mv "${OUT}/.work/${id}.exact" "${OUT}/exact/${id}.txt"
    echo "${id} ${len} $(cat "${OUT}/.work/${id}.time")" >> "${OUT}/exact_times.tsv"
    rm -f "${fa}" "${OUT}/.work/${id}.time"
}
export -f run_one
export OUT BIN
# long sequences first (slowest), then the Table 1 set
{ awk '{print $1, $2, "benchmark/data/bpRNA_1m_90_filtered_clean.fa"}' benchmark/data/bprna_long_1000_4381_list.txt
  awk '{print $1, $2, "benchmark/data/bprna_100_1000.fa"}' benchmark/data/bprna_100_1000_list.txt
} | xargs -P 16 -L 1 bash -c 'run_one "$@"' _
echo "End: $(date -u +%FT%TZ) w100=$(ls "${OUT}/w100" | wc -l) exact=$(ls "${OUT}/exact" | wc -l)"
