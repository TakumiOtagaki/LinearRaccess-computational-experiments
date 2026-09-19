#!/usr/bin/env python3
"""Rerun Raccess Exact at K=4 (access_len=1,3,7,10), matching the K=4 used
by the LinearRaccess/Raccess Window matched rerun (benchmark/results/
bprna_matched_runtime), so all three methods in Fig. 3 use the same K.
Raccess Exact still uses max_span=<sequence length> (no span restriction);
that is what makes it "Exact", not a K setting.

Usage:
    python3 run_raccess_exact_chunk_matched.py --batch BATCH.txt \
        --fasta bpRNA_1m_90.fasta --outdir OUT --raccess PATH --task-id N
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

ACCESS_LEN = "1,3,7,10"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_fasta_index(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name, parts = None, []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(parts)
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line)
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs


def parse_time_v(text: str) -> tuple[float | None, float | None, bool]:
    elapsed = None
    m = re.search(r"Elapsed \(wall clock\) time.*?(\d+):([\d.]+)$", text, re.MULTILINE)
    if m:
        elapsed = int(m.group(1)) * 60 + float(m.group(2))
    else:
        m = re.search(r"Elapsed \(wall clock\) time.*?(\d+):(\d+):([\d.]+)$", text, re.MULTILINE)
        if m:
            elapsed = int(m.group(1)) * 3600 + int(m.group(2)) * 60 + float(m.group(3))
    mem_kb = None
    m = re.search(r"Maximum resident set size.*?(\d+)", text)
    if m:
        mem_kb = int(m.group(1))
    ok = "Exit status: 0" in text
    return elapsed, mem_kb, ok


def run_timed(cmd: list[str], timeout: float) -> tuple[float | None, float | None, bool]:
    with tempfile.NamedTemporaryFile(mode="r", suffix=".time", delete=False) as tf:
        log_path = Path(tf.name)
    with log_path.open("w") as handle:
        subprocess.run(["/usr/bin/time", "-v", *cmd], stdout=subprocess.DEVNULL,
                        stderr=handle, timeout=timeout)
    text = log_path.read_text()
    log_path.unlink()
    return parse_time_v(text)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--raccess", type=Path, required=True)
    parser.add_argument("--task-id", type=int, required=True)
    parser.add_argument("--timeout", type=float, default=1800)
    args = parser.parse_args()

    batch_dir = args.outdir / f"batch_{args.task_id}"
    batch_dir.mkdir(parents=True, exist_ok=True)

    fasta_index = load_fasta_index(args.fasta)
    entries = []
    for line in args.batch.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        seq_id, length = line.split("\t")
        entries.append((seq_id, int(length)))

    rows = []
    for seq_id, length in entries:
        seq = fasta_index[seq_id]
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            seqfile = td / "seq.fa"
            seqfile.write_text(f">{seq_id}\n{seq}\n")
            out = td / "out.txt"
            cmd = [str(args.raccess), f"-seqfile={seqfile}", f"-access_len={ACCESS_LEN}",
                   f"-max_span={length}", f"-outfile={out}"]
            elapsed, mem_kb, ok = run_timed(cmd, args.timeout)
        rows.append(dict(seq_id=seq_id, length=length, elapsed_sec=elapsed,
                          peak_rss_mb=(mem_kb / 1024) if mem_kb else None, ok=ok))

    n_failed = sum(1 for r in rows if not r["ok"])
    manifest = dict(
        schema_version=1, task_id=args.task_id,
        finished_utc=datetime.now(timezone.utc).isoformat(),
        batch_sha256=digest(args.batch), raccess_sha256=digest(args.raccess),
        access_len=ACCESS_LEN, n_sequences=len(rows), n_failed=n_failed,
        records=rows,
    )
    (batch_dir / "run_summary.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Batch {args.task_id}: {len(rows)} sequences, {n_failed} failed", flush=True)
    return 0 if n_failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
