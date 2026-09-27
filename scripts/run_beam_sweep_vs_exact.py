#!/usr/bin/env python3
"""Beam sweep for LinearRaccess and LinearSampling against a common exact
Turner 2004 reference (docs/pi_feedback_20260925.md).

Reference (already cached by run_linracc_error_decomposition.py):
  LinearRaccess beam=0, c_hairpin=c_multi=len(seq), exact log-sum-exp.
Swept here, beam in BEAMS:
  LinearRaccess  beam=b, default caps (30), fast log-sum-exp (default)
  LinearSampling beam=b, k=K  (k large enough that beam error dominates)
beam=100 results are taken from run_exact_reference_si.py's cache.

Wall-clock times are recorded per run, but runs execute in parallel on one
node, so they are indicative only, not the paper's timing numbers.
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "benchmark/scripts"))
from run_raccess_compare import parse_linracc  # noqa: E402
from run_linearsampling_accuracy_vs_k import read_fasta_dict, select_sequences  # noqa: E402

LINRACC = REPO_ROOT / "vendor/LinearRaccess/build-raccess/LinRacc"
LINSAMP = REPO_ROOT / "vendor/LinearSampling/linearsampling"
COUNTER = REPO_ROOT / "benchmark/scripts/window_freq_counter"
FASTA = REPO_ROOT / "benchmark/data/bprna_100_1000.fa"
WORK_DIR = REPO_ROOT / "benchmark/results/.work/beam_sweep_vs_exact"
ACCESS_LENS = [1, 3, 7, 10]
BEAMS = [200, 300, 500]
K = 1_000_000


def run_linracc(name, seq, beam):
    with tempfile.TemporaryDirectory() as td:
        fa, out = Path(td) / "s.fa", Path(td) / "s.out"
        fa.write_text(f">{name}\n{seq}\n")
        t0 = time.perf_counter()
        subprocess.run([str(LINRACC), f"-seqfile={fa}", f"-outfile={out}",
                        f"-access_len={','.join(map(str, ACCESS_LENS))}", f"-beam={beam}",
                        "-energy=turner2004", "-probabilities"],
                       check=True, capture_output=True, text=True)
        sec = time.perf_counter() - t0
        parsed = parse_linracc(out)[name]
    return {str(w): [parsed[p][w] for p in range(len(seq) - w + 1)] for w in ACCESS_LENS}, sec


def run_linsamp(seq, beam):
    cmd = (f"{LINSAMP} -k {K} -b {beam} | tail -n +2 | "
           f"{COUNTER} {len(seq)} {','.join(map(str, ACCESS_LENS))}")
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, shell=True, input=seq + "\n", capture_output=True,
                          text=True, check=True)
    sec = time.perf_counter() - t0
    res = {str(w): [0.0] * (len(seq) - w + 1) for w in ACCESS_LENS}
    for line in proc.stdout.splitlines():
        w, pos, f = line.split()
        res[w][int(pos)] = float(f)
    return res, sec


def task(name, seq, tool, beam):
    cache = WORK_DIR / f"{name}__{tool}_b{beam}.json"
    if cache.exists():
        return
    if tool == "linracc":
        vals, sec = run_linracc(name, seq, beam)
    else:
        vals, sec = run_linsamp(seq, beam)
    cache.write_text(json.dumps({"sec": sec, "vals": vals}))


def main():
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    fasta = read_fasta_dict(FASTA)
    chosen = select_sequences(13)
    tasks = [(n, t, b) for n, _ in chosen for t in ("linracc", "linsamp") for b in BEAMS]
    tasks.sort(key=lambda x: -(len(fasta[x[0]]) * x[2] ** 2))  # heaviest first
    with ThreadPoolExecutor(max_workers=16) as ex:
        futs = [ex.submit(task, n, fasta[n], t, b) for n, t, b in tasks]
        for i, f in enumerate(as_completed(futs), 1):
            f.result()
            if i % 10 == 0:
                print(f"{i}/{len(tasks)}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
