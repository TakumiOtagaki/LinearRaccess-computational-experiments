#!/usr/bin/env python3
"""Supplementary data for the LinearSampling comparison, against an exact
Turner 2004 reference (docs/pi_feedback_20260925.md, "Revised plan").

For the same 52 bpRNA sequences / Table 1 length bins used elsewhere:
  reference : LinearRaccess beam=0 (exact), c_hairpin=c_multi=len(seq)
  compared  : LinearRaccess beam=100 (default caps 30),
              LinearSampling beam=100 at k in K_VALUES
Outputs, per length bin and window length w:
  - median/p95/max absolute error of each compared method vs. the reference
  - distribution of the exact accessibility values (quantiles, fraction below
    thresholds), and the implied LinearSampling sample count
    k ~ (1-p)/(p*eps^2) for relative standard error eps.

LinearRaccess and LinearSampling's exact values were checked to agree to
within 1e-3 (benchmark/results/worst_cell_diagnosis.json), so LinearRaccess
beam=0 serves as the common reference for both tools.

Sequences run in parallel (accuracy only; no timing is reported here).
Per-sequence results are cached under WORK_DIR so the job can be resumed.
"""
from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "benchmark/scripts"))
from run_raccess_compare import parse_linracc  # noqa: E402
from run_linearsampling_accuracy_vs_k import (  # noqa: E402
    BIN_LABELS, bin_for_length, read_fasta_dict, select_sequences,
)

LINRACC = REPO_ROOT / "vendor/LinearRaccess/build-raccess/LinRacc"
LINSAMP = REPO_ROOT / "vendor/LinearSampling/linearsampling"
COUNTER = REPO_ROOT / "benchmark/scripts/window_freq_counter"
FASTA = REPO_ROOT / "benchmark/data/bprna_100_1000.fa"
WORK_DIR = REPO_ROOT / "benchmark/results/.work/exact_reference_si"
ACCESS_LENS = [1, 3, 7, 10]
K_VALUES = [50_000, 1_000_000]
THRESHOLDS = [1e-3, 1e-5, 1e-7]
EPS_VALUES = [0.1, 0.01]


def linracc(name, seq, beam, caps, tmpdir):
    fa = tmpdir / f"{name}_b{beam}.fa"
    out = tmpdir / f"{name}_b{beam}.out"
    fa.write_text(f">{name}\n{seq}\n")
    cmd = [str(LINRACC), f"-seqfile={fa}", f"-outfile={out}",
           f"-access_len={','.join(map(str, ACCESS_LENS))}", f"-beam={beam}",
           "-energy=turner2004", "-probabilities"]
    if caps is not None:
        cmd += [f"-c_hairpin={caps}", f"-c_multi={caps}"]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    parsed = parse_linracc(out)[name]
    return {w: [parsed[pos][w] for pos in range(len(seq) - w + 1)] for w in ACCESS_LENS}


def linsamp(seq, k):
    cmd = (f"{LINSAMP} -k {k} -b 100 | tail -n +2 | "
           f"{COUNTER} {len(seq)} {','.join(map(str, ACCESS_LENS))}")
    proc = subprocess.run(cmd, shell=True, input=seq + "\n", capture_output=True,
                          text=True, check=True)
    res = {w: [0.0] * (len(seq) - w + 1) for w in ACCESS_LENS}
    for line in proc.stdout.splitlines():
        w, pos, f = line.split()
        res[int(w)][int(pos)] = float(f)
    return res


def process(name, seq):
    cache = WORK_DIR / f"{name}.json"
    if cache.exists():
        return json.loads(cache.read_text())
    with tempfile.TemporaryDirectory() as td:
        tmpdir = Path(td)
        result = {
            "name": name, "length": len(seq),
            "exact": linracc(name, seq, 0, len(seq), tmpdir),
            "linracc_b100": linracc(name, seq, 100, None, tmpdir),
        }
    for k in K_VALUES:
        result[f"linsamp_b100_k{k}"] = linsamp(seq, k)
    cache.write_text(json.dumps(result))
    return result


def quantile(sorted_vals, q):
    return sorted_vals[int(q * (len(sorted_vals) - 1))]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-per-bin", type=int, default=13)
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--out", type=Path, default=REPO_ROOT / "benchmark/results/exact_reference_si.json")
    args = parser.parse_args()

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    fasta = read_fasta_dict(FASTA)
    chosen = select_sequences(args.n_per_bin)
    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(process, name, fasta[name]): name for name, _ in chosen}
        for fut in as_completed(futs):
            r = fut.result()
            results.append(r)
            print(f"done {r['name']} len={r['length']} ({len(results)}/{len(chosen)})",
                  file=sys.stderr, flush=True)

    methods = ["linracc_b100"] + [f"linsamp_b100_k{k}" for k in K_VALUES]
    errors = {m: {b: {w: [] for w in ACCESS_LENS} for b in BIN_LABELS} for m in methods}
    exact_vals = {b: {w: [] for w in ACCESS_LENS} for b in BIN_LABELS}
    for r in results:
        b = BIN_LABELS[bin_for_length(r["length"])]
        for w in ACCESS_LENS:
            ex_w = r["exact"][str(w)] if str(w) in r["exact"] else r["exact"][w]
            exact_vals[b][w].extend(ex_w)
            for m in methods:
                m_w = r[m][str(w)] if str(w) in r[m] else r[m][w]
                errors[m][b][w].extend(abs(a - e) for a, e in zip(m_w, ex_w))

    summary = {"errors": {}, "exact_distribution": {}}
    for m in methods:
        summary["errors"][m] = {}
        for b in BIN_LABELS:
            summary["errors"][m][b] = {}
            for w in ACCESS_LENS:
                e = sorted(errors[m][b][w])
                summary["errors"][m][b][str(w)] = {
                    "n_cells": len(e), "median": statistics.median(e),
                    "p95": quantile(e, 0.95), "max": e[-1],
                }
    for b in BIN_LABELS:
        summary["exact_distribution"][b] = {}
        for w in ACCESS_LENS:
            v = sorted(exact_vals[b][w])
            d = {
                "n_cells": len(v),
                "p25": quantile(v, 0.25), "median": statistics.median(v), "p75": quantile(v, 0.75),
                "frac_below": {str(t): sum(1 for x in v if x < t) / len(v) for t in THRESHOLDS},
            }
            d["implied_k"] = {
                f"eps{eps}": {
                    q: ((1 - p) / (p * eps ** 2) if p > 0 else None)
                    for q, p in (("p25", d["p25"]), ("median", d["median"]))
                }
                for eps in EPS_VALUES
            }
            summary["exact_distribution"][b][str(w)] = d

    args.out.write_text(json.dumps({
        "n_sequences": len(results), "k_values": K_VALUES,
        "reference": "LinearRaccess beam=0, c_hairpin=c_multi=len(seq), Turner 2004",
        "summary": summary,
    }, indent=2))
    print(f"Wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
