#!/usr/bin/env python3
"""Diagnose the extreme LinearRaccess-vs-LinearSampling disagreements found
in linearsampling_accuracy_vs_k.json (max abs error ~0.99 at k=10M, in both
directions). For each flagged (sequence, w, pos) cell, compare four settings
against each other, all Turner 2004:

  LinearRaccess  beam=0   (exact), c_hairpin=c_multi=len(seq)
  LinearRaccess  beam=100,          c_hairpin=c_multi=len(seq)
  LinearSampling beam=0   (exact sampling), k=K_EXACT
  LinearSampling beam=100,                   k=K_BEAM

If the two beam=0 values agree, the disagreement is a beam-pruning effect
(and whichever beam=100 value is far from them is the one failing); if the
two beam=0 values disagree, it is a model/implementation difference or bug.
c_internal is 30 in both tools (hardcoded in LinearRaccess; SINGLE_MAX_LEN
in LinearSampling), so it is not varied.
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "benchmark/scripts"))
from run_raccess_compare import parse_linracc  # noqa: E402
from run_linearsampling_accuracy_vs_k import read_fasta_dict  # noqa: E402

LINRACC = REPO_ROOT / "vendor/LinearRaccess/build-raccess/LinRacc"
LINSAMP = REPO_ROOT / "vendor/LinearSampling/linearsampling"
COUNTER = REPO_ROOT / "benchmark/scripts/window_freq_counter"
FASTA = REPO_ROOT / "benchmark/data/bprna_100_1000.fa"
ACCESS_LENS = [1, 3, 7, 10]
K_EXACT = 1_000_000
K_BEAM = 10_000_000

CELLS = [
    ("bpRNA_CRW_2027", 1, 254),
    ("bpRNA_CRW_2030", 1, 14),
    ("bpRNA_CRW_1994", 3, 44),
]


def linracc(name, seq, beam, tmpdir):
    fa = tmpdir / f"{name}.fa"
    out = tmpdir / f"{name}_b{beam}.out"
    fa.write_text(f">{name}\n{seq}\n")
    cap = str(len(seq))
    subprocess.run([str(LINRACC), f"-seqfile={fa}", f"-outfile={out}",
                    f"-access_len={','.join(map(str, ACCESS_LENS))}", f"-beam={beam}",
                    "-energy=turner2004", f"-c_hairpin={cap}", f"-c_multi={cap}",
                    "-probabilities"], check=True, capture_output=True, text=True)
    return parse_linracc(out)[name]


def linsamp(seq, beam, k):
    cmd = (f"{LINSAMP} -k {k} -b {beam} | tail -n +2 | "
           f"{COUNTER} {len(seq)} {','.join(map(str, ACCESS_LENS))}")
    proc = subprocess.run(cmd, shell=True, input=seq + "\n", capture_output=True,
                          text=True, check=True)
    res = {}
    for line in proc.stdout.splitlines():
        w, pos, f = line.split()
        res[(int(w), int(pos))] = float(f)
    return res


def main():
    fasta = read_fasta_dict(FASTA)
    rows = []
    with tempfile.TemporaryDirectory() as td:
        tmpdir = Path(td)
        for name, w, pos in CELLS:
            seq = fasta[name]
            r = {"name": name, "length": len(seq), "w": w, "pos": pos}
            r["linracc_beam0"] = linracc(name, seq, 0, tmpdir)[pos][w]
            print(f"{name}: linracc beam0 done", file=sys.stderr, flush=True)
            r["linracc_beam100"] = linracc(name, seq, 100, tmpdir)[pos][w]
            r["linsamp_beam0"] = linsamp(seq, 0, K_EXACT)[(w, pos)]
            print(f"{name}: linsamp beam0 done", file=sys.stderr, flush=True)
            r["linsamp_beam100"] = linsamp(seq, 100, K_BEAM)[(w, pos)]
            print(json.dumps(r), flush=True)
            rows.append(r)
    out = REPO_ROOT / "benchmark/results/worst_cell_diagnosis.json"
    out.write_text(json.dumps({"k_exact": K_EXACT, "k_beam": K_BEAM, "rows": rows}, indent=2))


if __name__ == "__main__":
    main()
