#!/bin/bash
#$ -cwd
#$ -l cpu_16=1
#$ -l h_rt=12:00:00
#$ -N ls_long_calib
#$ -o logs/uge/linearsampling_comparison_long_calib_$JOB_ID.out
#$ -j y
#
# Cost calibration for the long-sequence (> beam) LinearSampling comparison:
# wall time and peak RSS of the exact reference (LinearRaccess beam=0, caps
# auto = N, exact log-sum-exp, Turner 2004, w in {1,3,7,10,35}) on one bpRNA
# sequence each of ~1,000 / 1,500 / 2,000 / 3,000 nt, run concurrently.
set -euo pipefail
PROJECT_ROOT=/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper
cd "${PROJECT_ROOT}"
W=benchmark/results/linearsampling_comparison_long/calib
BIN=benchmark/results/linearsampling_comparison_capb/bin/LinRacc
echo "Start: $(date -u +%FT%TZ) host=$(hostname)"
while read -r id len; do
  (
    awk -v id=">${id}" '/^>/ {keep = ($1 == id); if (keep) print $1; next} keep {printf "%s", $0} END {print ""}' \
        benchmark/data/bpRNA_1m_90_filtered_clean.fa > "${W}/${id}.fa"
    /usr/bin/time -v "${BIN}" -seqfile="${W}/${id}.fa" -outfile="${W}/${id}.out" \
        -access_len=1,3,7,10,35 -beam=0 -energy=turner2004 -c_multi=auto -c_hairpin=auto \
        -no-fast-logsumexp -probabilities > /dev/null 2> "${W}/${id}.time"
    echo "${id} ${len} $(grep -E 'Elapsed|Maximum resident' "${W}/${id}.time" | awk '{print $NF}' | tr '\n' ' ')"
  ) &
done < "${W}/list.tsv"
wait
echo "End: $(date -u +%FT%TZ)"
