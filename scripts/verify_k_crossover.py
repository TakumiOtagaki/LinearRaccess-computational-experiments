#!/usr/bin/env python3
"""Empirically verify (not extrapolate) whether LinearSampling's wall-clock
time exceeds LinearRaccess's at larger k, for a few representative lengths.
The earlier calibration (calibrate_linearsampling_k.py) only tested k up to
1,000,000, where LinearSampling was still faster than LinearRaccess at every
tested length; this checks larger k values directly.
"""
from __future__ import annotations

import json
import platform
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO_ROOT = Path("/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper")
LINRACC = REPO_ROOT / "vendor/LinearRaccess/build/LinRacc"
LINSAMP = REPO_ROOT / "vendor/LinearSampling/linearsampling"
FASTA = REPO_ROOT / "benchmark/data/bpRNA_1m_90_filtered_clean.fa"
BEAM = 100
K_VALUES = [2_000_000, 5_000_000, 10_000_000]
TARGET_LENGTHS = [400, 1000, 4000]
REPS = 2


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


def pick_representative(seqs, targets):
    chosen, used = [], set()
    for t in targets:
        best = min((s for s in seqs if s[0] not in used), key=lambda s: abs(len(s[1]) - t))
        chosen.append(best)
        used.add(best[0])
    return chosen


def time_linracc(name: str, seq: str, tmpdir: Path) -> float:
    seqfile = tmpdir / f"{name}.fa"
    outfile = tmpdir / f"{name}.out"
    seqfile.write_text(f">{name}\n{seq}\n")
    cmd = [str(LINRACC), f"-seqfile={seqfile}", f"-outfile={outfile}",
           "-access_len=1,3,7,10", f"-beam={BEAM}", "-energy=turner2004", "-probabilities"]
    t0 = time.perf_counter()
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return time.perf_counter() - t0


def time_linsamp(seq: str, k: int) -> float:
    t0 = time.perf_counter()
    subprocess.run([str(LINSAMP), "-k", str(k), "-b", str(BEAM)],
                   input=seq + "\n", capture_output=True, text=True, check=True)
    return time.perf_counter() - t0


def main() -> int:
    seqs = read_fasta_all(FASTA)
    chosen = pick_representative(seqs, TARGET_LENGTHS)
    rows = []
    with tempfile.TemporaryDirectory() as td:
        tmpdir = Path(td)
        for name, seq in chosen:
            length = len(seq)
            lr_times = [time_linracc(name, seq, tmpdir) for _ in range(REPS)]
            row = {"name": name, "length": length, "linracc_sec": lr_times, "linsamp_sec": {}}
            print(f"{name} len={length} linracc={lr_times}", file=sys.stderr, flush=True)
            for k in K_VALUES:
                ls_times = [time_linsamp(seq, k) for _ in range(REPS)]
                row["linsamp_sec"][str(k)] = ls_times
                print(f"  k={k} times={ls_times}", file=sys.stderr, flush=True)
            rows.append(row)
    out = REPO_ROOT / "benchmark/results/k_crossover_verification.json"
    out.write_text(json.dumps({"hostname": platform.node(), "beam": BEAM, "rows": rows}, indent=2))
    print(f"Wrote {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
