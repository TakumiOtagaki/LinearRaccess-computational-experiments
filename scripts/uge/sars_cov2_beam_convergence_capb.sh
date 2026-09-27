#!/bin/bash
#$ -cwd
#$ -l h_rt=12:00:00
#$ -l cpu_16=1
#$ -N sars_conv_capb
#$ -o logs/uge/sars_conv_capb_$JOB_ID.out
#$ -j y
#
# Section 3.4 SARS-CoV-2 beam-convergence figure, rerun with loop caps tied
# to the beam.  Same command as sars_cov2_beam_convergence_2026-09-16
# (-engine=raccess, Turner 1999, -access_len=1, -no-normalize,
# -print-ensemble-energy), caps auto, beams 5..200 as before plus 300 and 500.
# One run at a time; wall time and peak RSS recorded per beam.

set -euo pipefail
REPO_ROOT=/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper
cd "${REPO_ROOT}"
OUT=benchmark/results/sars_cov2_beam_convergence_capb
BIN="${OUT}/bin/LinRacc"
mkdir -p "${OUT}/bin"
cp vendor/LinearRaccess/build-raccess/LinRacc "${BIN}"
"${BIN}" --build-info > "${OUT}/bin/build_info.json"
sha256sum "${BIN}" benchmark/data/NC_045512.2.fa > "${OUT}/bin/sha256.txt"
echo "Start: $(date -u +%FT%TZ) host=$(hostname)"
for beam in 5 10 25 50 100 200 300 500; do
    /usr/bin/time -v "${BIN}" -seqfile=benchmark/data/NC_045512.2.fa \
        -outfile="${OUT}/out_beam${beam}.txt" -access_len=1 -engine=raccess \
        -beam="${beam}" -c_multi=auto -c_hairpin=auto -energy=turner1999 \
        -no-normalize -metadata -print-ensemble-energy \
        > /dev/null 2> "${OUT}/time_beam${beam}.log"
    echo "beam=${beam} $(grep -m1 ensemble_free_energy "${OUT}/out_beam${beam}.txt") $(grep -E 'Elapsed|Maximum resident' "${OUT}/time_beam${beam}.log" | tr -s ' ' | tr '\n' ' ')"
done
echo "End: $(date -u +%FT%TZ)"
