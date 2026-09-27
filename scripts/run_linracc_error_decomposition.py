#!/usr/bin/env python3
"""Decompose LinearRaccess's error (Turner 2004) on the same 52 bpRNA
sequences as run_exact_reference_si.py into beam, loop-cap, and fast
log-sum-exp contributions, by running four extra settings:

  A  beam=100, caps=len, fast LSE
  B  beam=0,   caps=30,  fast LSE
  C  beam=0,   caps=len, exact LSE   (numerically exact reference)
  D  beam=100, caps=30,  exact LSE

Already cached by run_exact_reference_si.py (fast LSE):
  beam=0/caps=len ("exact"), beam=100/caps=30 ("linracc_b100").
Per-sequence results are cached under WORK_DIR; the analysis step is
separate (--analyze) so it can be rerun without recomputation.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "benchmark/scripts"))
from run_raccess_compare import parse_linracc  # noqa: E402
from run_linearsampling_accuracy_vs_k import read_fasta_dict, select_sequences  # noqa: E402

LINRACC = REPO_ROOT / "vendor/LinearRaccess/build-raccess/LinRacc"
FASTA = REPO_ROOT / "benchmark/data/bprna_100_1000.fa"
WORK_DIR = REPO_ROOT / "benchmark/results/.work/linracc_error_decomposition"
ACCESS_LENS = [1, 3, 7, 10]
SETTINGS = {
    "A_b100_caplen_fast": (100, "len", True),
    "B_b0_cap30_fast": (0, "30", True),
    "C_b0_caplen_exactlse": (0, "len", False),
    "D_b100_cap30_exactlse": (100, "30", False),
}


def linracc(name, seq, beam, cap, fast, tmpdir, tag):
    fa = tmpdir / f"{name}_{tag}.fa"
    out = tmpdir / f"{name}_{tag}.out"
    fa.write_text(f">{name}\n{seq}\n")
    capv = str(len(seq)) if cap == "len" else cap
    cmd = [str(LINRACC), f"-seqfile={fa}", f"-outfile={out}",
           f"-access_len={','.join(map(str, ACCESS_LENS))}", f"-beam={beam}",
           "-energy=turner2004", f"-c_hairpin={capv}", f"-c_multi={capv}", "-probabilities"]
    if not fast:
        cmd.append("-no-fast-logsumexp")
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    parsed = parse_linracc(out)[name]
    return {str(w): [parsed[pos][w] for pos in range(len(seq) - w + 1)] for w in ACCESS_LENS}


def process(name, seq, tag):
    cache = WORK_DIR / f"{name}__{tag}.json"
    if cache.exists():
        return name, tag
    beam, cap, fast = SETTINGS[tag]
    with tempfile.TemporaryDirectory() as td:
        res = linracc(name, seq, beam, cap, fast, Path(td), tag)
    cache.write_text(json.dumps(res))
    return name, tag


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-per-bin", type=int, default=13)
    parser.add_argument("--jobs", type=int, default=16)
    args = parser.parse_args()
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    fasta = read_fasta_dict(FASTA)
    chosen = select_sequences(args.n_per_bin)
    tasks = [(name, tag) for name, _ in chosen for tag in SETTINGS]
    # longest first, so the slowest exact runs start early
    tasks.sort(key=lambda t: -len(fasta[t[0]]))
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(process, name, fasta[name], tag) for name, tag in tasks]
        for i, fut in enumerate(as_completed(futs), 1):
            name, tag = fut.result()
            print(f"done {name} {tag} ({i}/{len(tasks)})", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
