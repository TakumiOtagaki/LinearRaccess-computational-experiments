#!/usr/bin/env python3
"""Detailed (per-w, finer-quantile) version of the ad hoc cap-sensitivity
check: does raising LinearRaccess's c_hairpin/c_multi from 30 to the
sequence length change its own output, and if so how much and on which
window lengths? See docs/pi_feedback_20260925.md for the coarser first pass
that motivated this (negligible up to 2,000 nt, growing by 4,000 nt).
"""
from __future__ import annotations

import json
import platform
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "benchmark/scripts"))
from run_raccess_compare import parse_linracc  # noqa: E402

LINRACC = REPO_ROOT / "vendor/LinearRaccess/build-raccess/LinRacc"
FASTA = REPO_ROOT / "benchmark/data/bpRNA_1m_90_filtered_clean.fa"
ACCESS_LENS = [1, 3, 7, 10]
BEAM = 100
TARGET_LENGTHS = [150, 400, 700, 1000, 2000, 4000]


def read_fasta_all(path: Path):
    seqs, name, buf = [], None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs.append((name, "".join(buf)))
            name = line[1:].strip(); buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs.append((name, "".join(buf)))
    return seqs


def pick_representative(seqs, targets):
    chosen, used = [], set()
    for t in targets:
        best = min((s for s in seqs if s[0] not in used), key=lambda s: abs(len(s[1]) - t))
        chosen.append(best); used.add(best[0])
    return chosen


def run_linracc(name, seq, cap, tmpdir):
    seqfile = tmpdir / f"{name}_{cap}.fa"
    outfile = tmpdir / f"{name}_{cap}.out"
    seqfile.write_text(f">{name}\n{seq}\n")
    cmd = [str(LINRACC), f"-seqfile={seqfile}", f"-outfile={outfile}",
           f"-access_len={','.join(map(str, ACCESS_LENS))}", f"-beam={BEAM}",
           "-energy=turner2004", f"-c_hairpin={cap}", f"-c_multi={cap}", "-probabilities"]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return parse_linracc(outfile)[name]


def main():
    seqs = read_fasta_all(FASTA)
    chosen = pick_representative(seqs, TARGET_LENGTHS)
    rows = []
    with tempfile.TemporaryDirectory() as td:
        tmpdir = Path(td)
        for name, seq in chosen:
            length = len(seq)
            capped = run_linracc(name, seq, 30, tmpdir)
            uncapped = run_linracc(name, seq, length, tmpdir)
            per_w = {w: [] for w in ACCESS_LENS}
            for pos, vals in capped.items():
                for w, p in vals.items():
                    p2 = uncapped.get(pos, {}).get(w)
                    if p2 is not None:
                        per_w[w].append(abs(p - p2))
            row = {"name": name, "length": length, "per_w": {}}
            for w in ACCESS_LENS:
                diffs = sorted(per_w[w])
                n = len(diffs)
                if n == 0:
                    row["per_w"][str(w)] = None
                    continue
                row["per_w"][str(w)] = {
                    "n_cells": n,
                    "n_changed": sum(1 for d in diffs if d > 1e-9),
                    "median": statistics.median(diffs),
                    "p90": diffs[int(0.90 * (n - 1))],
                    "p99": diffs[int(0.99 * (n - 1))],
                    "p999": diffs[int(0.999 * (n - 1))],
                    "max": diffs[-1],
                }
                print(f"len={length:5d} w={w:2d} n={n:6d} median={row['per_w'][str(w)]['median']:.3e} "
                      f"p90={row['per_w'][str(w)]['p90']:.3e} p99={row['per_w'][str(w)]['p99']:.3e} "
                      f"max={row['per_w'][str(w)]['max']:.3e}", flush=True)
            rows.append(row)

    out = REPO_ROOT / "benchmark/results/linracc_cap_sensitivity.json"
    out.write_text(json.dumps({"hostname": platform.node(), "beam": BEAM, "rows": rows}, indent=2))
    print(f"Wrote {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
