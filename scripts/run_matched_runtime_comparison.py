#!/usr/bin/env python3
"""Matched-condition runtime comparison: LinearRaccess vs Raccess Window.

Fixes the three confounds codex exec flagged in the original Fig. 3 /
S3.3 "roughly 6-fold" comparison:
  1. Same access_len list (K=4: 1,3,7,10) for both tools, not K=4 vs K=1.
  2. Same timing boundary for both: /usr/bin/time -v around the binary
     invocation only, excluding Python-side I/O/validation/hashing.
  3. Strictly sequential (no thread pool / no concurrent array tasks),
     so no CPU contention confound.
"""
import json
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path("/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper")
LINRACC = ROOT / "vendor/LinearRaccess/build-raccess/LinRacc"
RACCESS = ROOT / "vendor/raccess/src/raccess/run_raccess"
SAMPLE_LIST = Path("/tmp/matched_bench/sample_list.tsv")
SAMPLE_FA = Path("/tmp/matched_bench/sample.fa")
OUT = Path("/tmp/matched_bench/results.jsonl")

ACCESS_LEN = "1,3,7,10"

def load_seqs():
    seqs = {}
    name = None
    parts = []
    for line in SAMPLE_FA.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(parts)
            name = line[1:].strip()
            parts = []
        else:
            parts.append(line.strip())
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs


def parse_time_v(stderr_text):
    elapsed = None
    m = re.search(r'Elapsed \(wall clock\) time.*?(\d+):([\d.]+)', stderr_text)
    if m:
        minutes, seconds = int(m.group(1)), float(m.group(2))
        elapsed = minutes * 60 + seconds
    else:
        m = re.search(r'Elapsed \(wall clock\) time.*?(\d+):(\d+):([\d.]+)', stderr_text)
        if m:
            hours, minutes, seconds = int(m.group(1)), int(m.group(2)), float(m.group(3))
            elapsed = hours * 3600 + minutes * 60 + seconds
    mem_kb = None
    m = re.search(r'Maximum resident set size.*?(\d+)', stderr_text)
    if m:
        mem_kb = int(m.group(1))
    exit_ok = "Exit status: 0" in stderr_text
    return elapsed, mem_kb, exit_ok


def run_timed(cmd):
    with tempfile.NamedTemporaryFile(mode="r", suffix=".time", delete=False) as tf:
        time_log = Path(tf.name)
    full_cmd = ["/usr/bin/time", "-v", *cmd]
    with time_log.open("w") as handle:
        subprocess.run(full_cmd, stdout=subprocess.DEVNULL, stderr=handle, timeout=600)
    text = time_log.read_text()
    time_log.unlink()
    return parse_time_v(text)


def main():
    seqs = load_seqs()
    order = [line.split("\t")[0] for line in SAMPLE_LIST.read_text().splitlines() if line.strip()]
    done = set()
    if OUT.exists():
        for line in OUT.read_text().splitlines():
            if line.strip():
                done.add(json.loads(line)["seq_id"])
    total = len(order)
    with OUT.open("a") as out:
        for i, seq_id in enumerate(order, 1):
            if seq_id in done:
                continue
            seq = seqs[seq_id]
            length = len(seq)
            with tempfile.TemporaryDirectory() as td:
                td = Path(td)
                seqfile = td / "seq.fa"
                seqfile.write_text(f">{seq_id}\n{seq}\n")

                lin_out = td / "lin.txt"
                lin_cmd = [str(LINRACC), f"-seqfile={seqfile}", f"-outfile={lin_out}",
                           f"-access_len={ACCESS_LEN}", "-engine=lincapr", "-energy=turner2004",
                           "-beam=100", "-c_multi=30", "-c_hairpin=30", "-probabilities", "-metadata"]
                lin_elapsed, lin_mem_kb, lin_ok = run_timed(lin_cmd)

                rac_out = td / "rac.txt"
                rac_cmd = [str(RACCESS), f"-seqfile={seqfile}", f"-access_len={ACCESS_LEN}",
                           "-max_span=100", f"-outfile={rac_out}"]
                rac_elapsed, rac_mem_kb, rac_ok = run_timed(rac_cmd)

            row = dict(seq_id=seq_id, length=length,
                       linracc_elapsed_sec=lin_elapsed, linracc_mem_kb=lin_mem_kb, linracc_ok=lin_ok,
                       raccess_elapsed_sec=rac_elapsed, raccess_mem_kb=rac_mem_kb, raccess_ok=rac_ok)
            out.write(json.dumps(row) + "\n")
            out.flush()
            print(f"[{i}/{total}] {seq_id} len={length} lin={lin_elapsed} rac={rac_elapsed}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
