#!/usr/bin/env python3
"""Run one bpRNA chunk through both LinearRaccess and Raccess Window at
matched conditions (same K=4 access_len list, executable-only timing via
/usr/bin/time -v, strictly sequential), replacing the original Fig. 3/S3.3
production data (benchmark/results/bprna_runs, raccess_w100), which used
K=4 vs K=1 and different timing boundaries (see
benchmark/results/matched_runtime_comparison_20260919 for the smaller
controlled check that first identified this).

Usage:
    python3 run_bprna_chunk_matched.py --chunk CHUNK.fa --outdir OUT \
        --linracc PATH --raccess PATH --task-id N
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

ACCESS_LEN = "1,3,7,10"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records, name, parts = [], None, []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(parts)))
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line)
    if name is not None:
        records.append((name, "".join(parts)))
    return records


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
    parser.add_argument("--chunk", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--linracc", type=Path, required=True)
    parser.add_argument("--raccess", type=Path, required=True)
    parser.add_argument("--task-id", type=int, required=True)
    parser.add_argument("--timeout", type=float, default=1800)
    parser.add_argument("--access-len", default=ACCESS_LEN)
    parser.add_argument("--beam", type=int, default=100)
    parser.add_argument("--caps", default="30",
                        help="LinearRaccess c_multi/c_hairpin: an integer or 'auto' (= beam)")
    args = parser.parse_args()

    batch = args.outdir / f"batch_{args.task_id}"
    batch.mkdir(parents=True, exist_ok=True)

    linracc_info = json.loads(subprocess.run(
        [str(args.linracc), "--build-info"], capture_output=True, text=True, check=True).stdout)

    records = parse_fasta(args.chunk)
    rows = []
    for name, seq in records:
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            seqfile = td / "seq.fa"
            seqfile.write_text(f">{name}\n{seq}\n")

            lin_out = td / "lin.txt"
            lin_cmd = [str(args.linracc), f"-seqfile={seqfile}", f"-outfile={lin_out}",
                       f"-access_len={args.access_len}", "-engine=lincapr", "-energy=turner2004",
                       f"-beam={args.beam}", f"-c_multi={args.caps}", f"-c_hairpin={args.caps}",
                       "-probabilities", "-metadata"]
            lin_elapsed, lin_mem_kb, lin_ok = run_timed(lin_cmd, args.timeout)

            rac_out = td / "rac.txt"
            rac_cmd = [str(args.raccess), f"-seqfile={seqfile}", f"-access_len={args.access_len}",
                       "-max_span=100", f"-outfile={rac_out}"]
            rac_elapsed, rac_mem_kb, rac_ok = run_timed(rac_cmd, args.timeout)

        rows.append(dict(
            seq_id=name, length=len(seq),
            linracc_elapsed_sec=lin_elapsed,
            linracc_peak_rss_mb=(lin_mem_kb / 1024) if lin_mem_kb else None,
            linracc_ok=lin_ok,
            raccess_elapsed_sec=rac_elapsed,
            raccess_peak_rss_mb=(rac_mem_kb / 1024) if rac_mem_kb else None,
            raccess_ok=rac_ok,
        ))

    n_failed = sum(1 for r in rows if not r["linracc_ok"] or not r["raccess_ok"])
    manifest = dict(
        schema_version=1, task_id=args.task_id,
        finished_utc=datetime.now(timezone.utc).isoformat(),
        chunk_sha256=digest(args.chunk),
        linracc_sha256=digest(args.linracc), linracc_build_info=linracc_info,
        raccess_sha256=digest(args.raccess),
        n_sequences=len(rows), n_failed=n_failed,
        access_len=args.access_len,
        settings=dict(linearraccess=dict(beam=args.beam, c_multi=args.caps, c_hairpin=args.caps,
                                          engine="lincapr", energy="turner2004"),
                      raccess_window=dict(max_span=100)),
        records=rows,
    )
    (batch / "run_summary.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Batch {args.task_id}: {len(rows)} sequences, {n_failed} failed", flush=True)
    return 0 if n_failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
