#!/usr/bin/env python3
"""LinearSampling (k=50,000, own default beam=100, own compiled-in Turner
2004 parameters) vs. LinearRaccess (beam=100, energy=turner2004, own
defaults).

k=50,000 matches LinearSampling's own paper (Zhang et al. 2023)'s sample
size for its reported experiments on sequences of comparable length
(<=4,000 nt) and its SARS-CoV-2 application, rather than an arbitrary or
compute-budget-matched choice -- author decision 2026-09-25, after finding
that matching wall-clock time to LinearRaccess would require k on the order
of 1.5-4.5 million (benchmark/scripts/calibrate_linearsampling_k.py,
benchmark/results/linearsampling_k_calibration.json), and after noting that
a beam=100-vs-beam=100 point-estimate comparison can't cleanly isolate
sampling noise from each tool's own beam-approximation error anyway (that
would need an exact/beam=0 reference on both sides, which is a separate,
heavier experiment). This run's actual per-sequence wall-clock time is
still recorded for context, not as a claimed fairness property.

Loop caps: LinearRaccess's c_hairpin/c_multi are raised to a fixed 900 (not
the sequence length -- capping at N would make C scale with N, turning the
O(N(b^2+Kb(C+c_internal^2))) complexity quadratic and breaking the linear-
time property the whole paper is about; 900 is a fixed constant, same order
of magnitude as c_internal^2=30^2=900, so linearity in N is untouched) to
approximate LinearSampling's behavior, which does not actually enforce a
hairpin/multiloop length cap despite defining HAIRPIN_MAX_LEN/
MULTI_MAX_LEN=30 (those constants are unused in its .cpp sources; only the
internal/bulge cap, SINGLE_MAX_LEN=30, is enforced). LinearRaccess's own
internal/bulge cap (c_internal) has no CLI flag and stays at its hardcoded
30, already matching SINGLE_MAX_LEN.

A direct empirical check (benchmark/scripts (ad hoc, not committed) on
bpRNA sequences at 150-4000 nt) found the cap 30-vs-900 difference in
LinearRaccess's own output is negligible (median 1e-9 to 1e-8, max ~1e-5
to 9e-5) up to 2,000 nt, but grows to a max of 5.8e-3 by 4,000 nt --
comparable to Table 1's own error scale at that point. Within this
script's actual tested range (100-1,000 nt, matching Table 1's bins), the
effect stays negligible (max ~5.7e-5 at 1,000 nt), so the k=50,000 result
already computed with c_hairpin=c_multi=30 (Table 1's own convention)
did not need to be rerun. This 900 setting is kept here as the correct
default if this comparison is ever extended to longer sequences.

Mirrors Table 1's methodology: same bpRNA sequences (subset of
benchmark/data/bprna_100_1000_list.txt), same window lengths
w in {1,3,7,10}, same four length bins, median/p95/max absolute error
per bin. Unlike Table 1, the reference here is LinearRaccess itself (not
Raccess Exact), since LinearSampling's compiled-in energy parameters match
LinearRaccess's own Turner 2004 default, not Raccess's Turner 1999-only
support -- Raccess/Raccess Exact cannot run this comparison's energy model
at all.

Usage:
    python3 benchmark/scripts/run_linearsampling_vs_linearraccess.py \
        [--n-per-bin 13] [--out benchmark/results/linearsampling_accuracy.json]
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
K_SAMPLES = 50_000
BINS = [(100, 300), (300, 500), (500, 700), (700, 1001)]
BIN_LABELS = ["100-300", "300-500", "500-700", "700-1000"]


def bin_for_length(length: int) -> int | None:
    for idx, (lo, hi) in enumerate(BINS):
        if lo <= length < hi:
            return idx
    return None


def read_fasta_dict(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    buf: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            name = line[1:].strip()
            buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def select_sequences(n_per_bin: int) -> list[tuple[str, int]]:
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


def run_linracc(name: str, seq: str, tmpdir: Path) -> dict[int, dict[int, float]]:
    seqfile = tmpdir / f"{name}.fa"
    outfile = tmpdir / f"{name}.out"
    seqfile.write_text(f">{name}\n{seq}\n")
    # c_hairpin/c_multi raised to a fixed 900 (see module docstring) to
    # approximate LinearSampling's actual uncapped hairpin/multiloop behavior,
    # without making C scale with N (which would break linear-time scaling).
    cap = "900"
    cmd = [
        str(LINRACC), f"-seqfile={seqfile}", f"-outfile={outfile}",
        f"-access_len={','.join(map(str, ACCESS_LENS))}", f"-beam={BEAM}",
        "-energy=turner2004", f"-c_hairpin={cap}", f"-c_multi={cap}", "-probabilities",
    ]
    t0 = time.perf_counter()
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    parsed = parse_linracc(outfile)[name]
    result: dict[int, dict[int, float]] = {w: {} for w in ACCESS_LENS}
    for pos, values in parsed.items():
        for length, prob in values.items():
            if length in result:
                result[length][pos] = prob
    return result, elapsed


WINDOW_FREQ_COUNTER = REPO_ROOT / "benchmark/scripts/window_freq_counter"


def run_linearsampling_and_accumulate(seq: str, w_list: list[int], length: int):
    """Pipe LinearSampling directly into window_freq_counter (a small
    compiled C++ helper, benchmark/scripts/window_freq_counter.cpp), which
    does the per-sample maximal-unpaired-run scan and difference-array
    accumulation (same idea as this paper's own Algorithm 2). Three
    pure-Python versions of this were tried and measured on a ~480-nt,
    k=5,000,000 test case before settling on this: a naive per-window slice
    check (too slow to finish inside a 3-hour walltime), a per-character
    while-loop run-scan (264.6s), and an re.finditer-based run-scan
    (187.2s) -- all far slower than LinearSampling's own raw execution time
    (30.6s). This C++ pipeline takes 27.0s, i.e. Python-level per-sample
    work was the actual bottleneck at this k, not node speed or LinearSampling
    itself. Only orchestration (launching the pipeline, parsing its small
    (length-w+1)-per-w text output) stays in Python."""
    valid_ws = [w for w in w_list if length - w + 1 > 0]
    t0 = time.perf_counter()
    linsamp_proc = subprocess.Popen(
        [str(LINSAMP), "-k", str(K_SAMPLES), "-b", str(BEAM)],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE,
    )
    tail_proc = subprocess.Popen(["tail", "-n", "+2"], stdin=linsamp_proc.stdout, stdout=subprocess.PIPE)
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

    result = {w: [0.0] * (length - w + 1) for w in valid_ws}
    for line in out.splitlines():
        w_s, pos_s, freq_s = line.split()
        result[int(w_s)][int(pos_s)] = float(freq_s)
    return result, elapsed, n_samples


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-per-bin", type=int, default=13)
    parser.add_argument("--out", type=Path, default=REPO_ROOT / "benchmark/results/linearsampling_accuracy.json")
    args = parser.parse_args()

    fasta = read_fasta_dict(FASTA)
    chosen = select_sequences(args.n_per_bin)
    print(f"Selected {len(chosen)} sequences ({args.n_per_bin} per bin)", file=sys.stderr)

    per_bin_errors: dict[int, dict[int, list[float]]] = {
        b: {w: [] for w in ACCESS_LENS} for b in range(len(BINS))
    }
    linracc_times: list[float] = []
    linsamp_times: list[float] = []
    per_sequence: list[dict] = []

    with tempfile.TemporaryDirectory() as td:
        tmpdir = Path(td)
        for i, (name, length) in enumerate(chosen):
            seq = fasta[name]
            b = bin_for_length(length)
            print(f"[{i+1}/{len(chosen)}] {name} len={length} bin={BIN_LABELS[b]}", file=sys.stderr, flush=True)
            linracc_probs, lr_elapsed = run_linracc(name, seq, tmpdir)
            freqs_by_w, ls_elapsed, n_samples = run_linearsampling_and_accumulate(seq, ACCESS_LENS, length)
            linracc_times.append(lr_elapsed)
            linsamp_times.append(ls_elapsed)
            print(f"  linracc={lr_elapsed:.2f}s linsamp={ls_elapsed:.2f}s (n_samples={n_samples})",
                  file=sys.stderr, flush=True)
            seq_errors = {w: [] for w in ACCESS_LENS}
            for w, freqs in freqs_by_w.items():
                for pos, p_hat in enumerate(freqs):
                    p_lin = linracc_probs[w].get(pos)
                    if p_lin is None:
                        continue
                    err = abs(p_hat - p_lin)
                    per_bin_errors[b][w].append(err)
                    seq_errors[w].append(err)
            per_sequence.append({
                "name": name, "length": length, "bin": BIN_LABELS[b],
                "linracc_sec": lr_elapsed, "linsampling_sec": ls_elapsed,
                "n_samples": n_samples,
                "median_abs_err": {w: (statistics.median(e) if e else None) for w, e in seq_errors.items()},
            })

    summary = {"bins": {}}
    for b, label in enumerate(BIN_LABELS):
        summary["bins"][label] = {}
        for w in ACCESS_LENS:
            errs = per_bin_errors[b][w]
            if not errs:
                summary["bins"][label][str(w)] = None
                continue
            errs_sorted = sorted(errs)
            summary["bins"][label][str(w)] = {
                "n_cells": len(errs),
                "median": statistics.median(errs_sorted),
                "p95": errs_sorted[int(0.95 * (len(errs_sorted) - 1))],
                "max": errs_sorted[-1],
            }

    result = {
        "hostname": platform.node(),
        "beam": BEAM,
        "k_samples": K_SAMPLES,
        "n_per_bin": args.n_per_bin,
        "n_sequences": len(chosen),
        "linracc_total_sec": sum(linracc_times),
        "linsampling_total_sec": sum(linsamp_times),
        "linracc_mean_sec": statistics.mean(linracc_times),
        "linsampling_mean_sec": statistics.mean(linsamp_times),
        "summary": summary,
        "per_sequence": per_sequence,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2))
    print(f"Wrote {args.out}", file=sys.stderr)
    print(f"Total: linracc={sum(linracc_times):.1f}s linsampling={sum(linsamp_times):.1f}s", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
