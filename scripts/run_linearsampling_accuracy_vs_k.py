#!/usr/bin/env python3
"""LinearSampling accuracy (and wall-clock time) vs. LinearRaccess, at
several sample counts k, on the same bpRNA sequences used by
run_linearsampling_vs_linearraccess.py (k=50,000 only). This checks whether
accuracy has converged/plateaued by k=10,000,000 -- and whether that
convergence happens before or after LinearSampling's own time exceeds
LinearRaccess's (already established separately: the crossover happens
around k=2,000,000-5,000,000; see verify_k_crossover.py /
benchmark/results/k_crossover_verification.json).

Uses LinearRaccess with c_hairpin=c_multi=30 (Table 1's own convention;
the loop-cap mismatch with LinearSampling was checked separately and found
negligible in this 100-1,000 nt range -- see
benchmark/scripts/check_linracc_cap_sensitivity.py).

Usage:
    python3 benchmark/scripts/run_linearsampling_accuracy_vs_k.py \
        [--n-per-bin 13] [--k-values 50000,1000000,5000000,10000000] \
        [--out benchmark/results/linearsampling_accuracy_vs_k.json]
"""
from __future__ import annotations

import argparse
import json
import platform
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "benchmark/scripts"))
from run_raccess_compare import parse_linracc  # noqa: E402

LINRACC = REPO_ROOT / "vendor/LinearRaccess/build-raccess/LinRacc"
LINSAMP = REPO_ROOT / "vendor/LinearSampling/linearsampling"
FASTA = REPO_ROOT / "benchmark/data/bprna_100_1000.fa"
SEQ_LIST = REPO_ROOT / "benchmark/data/bprna_100_1000_list.txt"
ACCESS_LENS = [1, 3, 7, 10]
BEAM = 100
BINS = [(100, 300), (300, 500), (500, 700), (700, 1001)]
BIN_LABELS = ["100-300", "300-500", "500-700", "700-1000"]


def bin_for_length(length: int):
    for idx, (lo, hi) in enumerate(BINS):
        if lo <= length < hi:
            return idx
    return None


def read_fasta_dict(path: Path) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            name = line[1:].strip(); buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def select_sequences(n_per_bin: int):
    per_bin: list[list[tuple[str, int]]] = [[] for _ in BINS]
    for line in SEQ_LIST.read_text().splitlines():
        name, length_s = line.split("\t")
        length = int(length_s)
        b = bin_for_length(length)
        if b is not None and len(per_bin[b]) < n_per_bin:
            per_bin[b].append((name, length))
    chosen = []
    for b in per_bin:
        chosen.extend(b)
    return chosen


def run_linracc(name, seq, tmpdir):
    seqfile = tmpdir / f"{name}.fa"
    outfile = tmpdir / f"{name}.out"
    seqfile.write_text(f">{name}\n{seq}\n")
    cmd = [str(LINRACC), f"-seqfile={seqfile}", f"-outfile={outfile}",
           f"-access_len={','.join(map(str, ACCESS_LENS))}", f"-beam={BEAM}",
           "-energy=turner2004", "-probabilities"]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    parsed = parse_linracc(outfile)[name]
    result = {w: {} for w in ACCESS_LENS}
    for pos, values in parsed.items():
        for length, prob in values.items():
            if length in result:
                result[length][pos] = prob
    return result


WINDOW_FREQ_COUNTER = REPO_ROOT / "benchmark/scripts/window_freq_counter"


def run_and_accumulate(seq, k, w_list, length):
    """Run LinearSampling piped directly into window_freq_counter (a small
    compiled C++ helper, benchmark/scripts/window_freq_counter.cpp) which
    does the per-sample maximal-unpaired-run scan and difference-array
    accumulation. Three pure-Python versions of this loop were tried and
    measured, in order, on the same ~480-nt/k=5,000,000 test case:
      - naive per-window slice check:        far too slow to finish (hit a
        3-hour walltime kill on the real job)
      - per-character while-loop run-scan:   264.6s
      - re.finditer-based run-scan:          187.2s
      - this C++ pipeline:                    27.0s (~= LinearSampling's own
        raw execution time, 30.6s with no post-processing at all)
    Python-level per-sample/per-character work is simply too slow at
    k=1e6-1e7 scale; only orchestration (launching the pipeline, parsing the
    resulting small (length-w+1)-per-w output) stays in Python."""
    valid_ws = [w for w in w_list if length - w + 1 > 0]
    t0 = time.perf_counter()
    linsamp_proc = subprocess.Popen(
        [str(LINSAMP), "-k", str(k), "-b", str(BEAM)],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE,
    )
    tail_proc = subprocess.Popen(
        ["tail", "-n", "+2"], stdin=linsamp_proc.stdout, stdout=subprocess.PIPE,
    )
    linsamp_proc.stdout.close()
    counter_proc = subprocess.Popen(
        [str(WINDOW_FREQ_COUNTER), str(length), ",".join(map(str, valid_ws))],
        stdin=tail_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    tail_proc.stdout.close()
    linsamp_proc.stdin.write((seq + "\n").encode())
    linsamp_proc.stdin.close()
    out, err = counter_proc.communicate()
    elapsed = time.perf_counter() - t0
    for p in (linsamp_proc, tail_proc, counter_proc):
        p.wait()
        if p.returncode != 0:
            raise subprocess.CalledProcessError(p.returncode, p.args)
    n_samples = int(err.strip().rsplit("=", 1)[1])

    # window_freq_counter already prefix-sums and divides by n_samples; each
    # output line is "w pos freq" directly.
    result = {w: [0.0] * (length - w + 1) for w in valid_ws}
    for line in out.splitlines():
        w_s, pos_s, freq_s = line.split()
        result[int(w_s)][int(pos_s)] = float(freq_s)
    return result, elapsed, n_samples


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-per-bin", type=int, default=13)
    parser.add_argument("--k-values", type=str, default="50000,1000000,5000000,10000000")
    parser.add_argument("--out", type=Path, default=REPO_ROOT / "benchmark/results/linearsampling_accuracy_vs_k.json")
    args = parser.parse_args()
    k_values = [int(x) for x in args.k_values.split(",")]

    fasta = read_fasta_dict(FASTA)
    chosen = select_sequences(args.n_per_bin)
    print(f"Selected {len(chosen)} sequences; k_values={k_values}", file=sys.stderr)

    per_k_bin_errors = {k: {b: {w: [] for w in ACCESS_LENS} for b in range(len(BINS))} for k in k_values}
    per_k_total_time = {k: 0.0 for k in k_values}
    linracc_total_time = 0.0

    with tempfile.TemporaryDirectory() as td:
        tmpdir = Path(td)
        for i, (name, length) in enumerate(chosen):
            seq = fasta[name]
            b = bin_for_length(length)
            t0 = time.perf_counter()
            linracc_probs = run_linracc(name, seq, tmpdir)
            linracc_total_time += time.perf_counter() - t0
            print(f"[{i+1}/{len(chosen)}] {name} len={length} bin={BIN_LABELS[b]}", file=sys.stderr, flush=True)
            valid_ws = [w for w in ACCESS_LENS if length - w + 1 > 0]
            for k in k_values:
                freqs_by_w, elapsed, n_samples = run_and_accumulate(seq, k, valid_ws, length)
                per_k_total_time[k] += elapsed
                print(f"  k={k}: {elapsed:.2f}s (n_samples={n_samples})", file=sys.stderr, flush=True)
                for w, freqs in freqs_by_w.items():
                    for pos, p_hat in enumerate(freqs):
                        p_lin = linracc_probs[w].get(pos)
                        if p_lin is None:
                            continue
                        per_k_bin_errors[k][b][w].append(abs(p_hat - p_lin))

    summary = {}
    for k in k_values:
        summary[str(k)] = {"total_linsampling_sec": per_k_total_time[k], "bins": {}}
        for b, label in enumerate(BIN_LABELS):
            summary[str(k)]["bins"][label] = {}
            for w in ACCESS_LENS:
                errs = sorted(per_k_bin_errors[k][b][w])
                if not errs:
                    summary[str(k)]["bins"][label][str(w)] = None
                    continue
                n = len(errs)
                summary[str(k)]["bins"][label][str(w)] = {
                    "n_cells": n,
                    "median": statistics.median(errs),
                    "p95": errs[int(0.95 * (n - 1))],
                    "p99": errs[int(0.99 * (n - 1))],
                    "max": errs[-1],
                }

    result = {
        "hostname": platform.node(), "beam": BEAM, "n_per_bin": args.n_per_bin,
        "n_sequences": len(chosen), "k_values": k_values,
        "linracc_total_sec": linracc_total_time,
        "summary": summary,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2))
    print(f"Wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
