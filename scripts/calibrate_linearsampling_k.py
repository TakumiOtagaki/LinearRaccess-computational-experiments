#!/usr/bin/env python3
"""Calibrate LinearSampling's sample number k so that its wall-clock time
roughly matches LinearRaccess's, at each of several representative sequence
lengths. Both tools run single-threaded, own default parameters otherwise
(LinearRaccess: beam=100, energy=turner2004; LinearSampling: beam=100).

Context: LinearSampling has no built-in parallelism (no thread/job CLI flag)
and reuses one beam-pruned parse forest across all k samples, so its runtime
is nearly flat in k -- an arbitrary small k (e.g. its own default, 10) makes
it look far noisier than a fair, equal-compute-budget comparison would. This
script measures wall-clock time as a function of k so a later accuracy run
can use a per-length k that spends roughly the same wall-clock time as
LinearRaccess, rather than an arbitrary sample count.

Not run on a login/interactive node: earlier interactive-session timings for
LinearRaccess (~180s at 4,000 nt) did not match the frozen benchmark bundle's
reported rate (~0.0059 s/nt, i.e. ~24s at 4,000 nt) -- about 7x slower,
consistent with node contention distorting absolute times. Submit this via
benchmark/scripts/uge/calibrate_linearsampling_k.sh on a dedicated node
instead.

Usage:
    python3 benchmark/scripts/calibrate_linearsampling_k.py [--out PATH]
"""
from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
LINRACC = REPO_ROOT / "vendor/LinearRaccess/build/LinRacc"
LINSAMP = REPO_ROOT / "vendor/LinearSampling/linearsampling"
FASTA = REPO_ROOT / "benchmark/data/bpRNA_1m_90_filtered_clean.fa"
BEAM = 100
K_VALUES = [10, 1000, 50_000, 200_000, 1_000_000]
TARGET_LENGTHS = [150, 400, 700, 1000, 2000, 4000]
REPS = 3


def read_fasta_all(path: Path) -> list[tuple[str, str]]:
    seqs = []
    name = None
    buf: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs.append((name, "".join(buf)))
            name = line[1:].strip()
            buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs.append((name, "".join(buf)))
    return seqs


def pick_representative(seqs: list[tuple[str, str]], targets: list[int]) -> list[tuple[str, str]]:
    chosen = []
    used = set()
    for t in targets:
        best = min((s for s in seqs if s[0] not in used), key=lambda s: abs(len(s[1]) - t))
        chosen.append(best)
        used.add(best[0])
    return chosen


def time_linracc(name: str, seq: str, tmpdir: Path) -> float:
    seqfile = tmpdir / f"{name}.fa"
    outfile = tmpdir / f"{name}.out"
    seqfile.write_text(f">{name}\n{seq}\n")
    cmd = [
        str(LINRACC), f"-seqfile={seqfile}", f"-outfile={outfile}",
        "-access_len=1,3,7,10", f"-beam={BEAM}", "-energy=turner2004", "-probabilities",
    ]
    t0 = time.perf_counter()
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return time.perf_counter() - t0


def time_linsamp(seq: str, k: int) -> float:
    t0 = time.perf_counter()
    subprocess.run(
        [str(LINSAMP), "-k", str(k), "-b", str(BEAM)],
        input=seq + "\n", capture_output=True, text=True, check=True,
    )
    return time.perf_counter() - t0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=REPO_ROOT / "benchmark/results/linearsampling_k_calibration.json")
    args = parser.parse_args()

    seqs = read_fasta_all(FASTA)
    chosen = pick_representative(seqs, TARGET_LENGTHS)
    rows = []
    with tempfile.TemporaryDirectory() as td:
        tmpdir = Path(td)
        for name, seq in chosen:
            length = len(seq)
            lr_times = [time_linracc(name, seq, tmpdir) for _ in range(REPS)]
            row = {
                "name": name,
                "length": length,
                "linracc_beam100_turner2004_sec": lr_times,
                "linearsampling_beam100_sec": {},
            }
            for k in K_VALUES:
                ls_times = [time_linsamp(seq, k) for _ in range(REPS)]
                row["linearsampling_beam100_sec"][str(k)] = ls_times
                print(f"{name} len={length} k={k} times={ls_times}", file=sys.stderr, flush=True)
            rows.append(row)
            print(f"{name} len={length} linracc_times={lr_times}", file=sys.stderr, flush=True)

    result = {
        "hostname": platform.node(),
        "beam": BEAM,
        "reps": REPS,
        "k_values": K_VALUES,
        "target_lengths": TARGET_LENGTHS,
        "rows": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2))
    print(f"Wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
