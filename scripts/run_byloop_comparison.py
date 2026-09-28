#!/usr/bin/env python3
"""Loop-type-resolved window accessibility: LinearRaccess (-byloop) vs
LinearSampling (samples classified by window_freq_counter_byloop).

For each window [i, i+w-1] and loop type c in {E, H, B, I, M} (exterior,
hairpin, bulge, internal, multiloop), the probability that the window is
unpaired inside a maximal unpaired run of loop type c (Methods Eq. for
P_u^c).  LinearRaccess reports these as a by-product of the total; for
LinearSampling each sampled structure is parsed to assign its unpaired runs
a loop type.

Sets and energy model as run_linearsampling_comparison_capb.py (Turner 2004):
--set bprna52 (13 per Table 1 bin) or --set long (20 sequences, 1,000-3,541 nt).
--phase accuracy (parallel): exact (LinearRaccess beam=0, exact LSE),
  linracc_b{200,500}, linsamp_b{200,500}_k1000000.
--phase timing (sequential, /usr/bin/time -v around the whole command):
  the same configurations except exact.
Per-(sequence, configuration) caches under --work-dir; the frozen
LinearRaccess binary is <work-dir>/bin/LinRacc.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "benchmark/scripts"))
import run_linearsampling_comparison_capb as base  # noqa: E402

LINSAMP = REPO_ROOT / "vendor/LinearSampling/linearsampling"
COUNTER = REPO_ROOT / "benchmark/scripts/window_freq_counter_byloop"
ACCESS_LENS = [1, 3, 7, 10, 35]
TYPES = "EHBIM"
BEAMS = [200, 500]
K = 1_000_000


def configs(phase: str) -> list[str]:
    cfg = [f"linracc_b{b}" for b in BEAMS] + [f"linsamp_b{b}_k{K}" for b in BEAMS]
    return (["exact"] + cfg) if phase == "accuracy" else cfg


def parse_linracc_byloop(path: Path, n: int) -> dict:
    vals = {t: {str(w): [0.0] * (n - w + 1) for w in ACCESS_LENS} for t in TYPES}
    for line in path.read_text().splitlines():
        if not line or line[0] in "#>" or "byloop:" not in line:
            continue
        pos = int(line.split("\t", 1)[0])
        for block in line.split("byloop:", 1)[1].strip().strip(";").split(";"):
            w, rest = block.split(",", 1)
            for kv in rest.split(","):
                t, v = kv.split("=")
                vals[t][w][pos] = float(v)
    return vals


def parse_counter(stdout: str, n: int) -> dict:
    vals = {t: {str(w): [0.0] * (n - w + 1) for w in ACCESS_LENS} for t in TYPES}
    for line in stdout.splitlines():
        t, w, pos, f = line.split()
        vals[t][w][int(pos)] = float(f)
    return vals


def run_one(work: Path, linracc: Path, phase: str, name: str, seq: str, cfg: str) -> None:
    cache = work / phase / f"{name}__{cfg}.json"
    if cache.exists():
        return
    n = len(seq)
    lens = ",".join(map(str, ACCESS_LENS))
    with tempfile.TemporaryDirectory(dir=work / "tmp") as td:
        fa, out = Path(td) / "s.fa", Path(td) / "s.out"
        fa.write_text(f">{name}\n{seq}\n")
        if cfg.startswith("linsamp"):
            beam, k = map(int, re.fullmatch(r"linsamp_b(\d+)_k(\d+)", cfg).groups())
            shell = f"{LINSAMP} -k {k} -b {beam} | tail -n +2 | {COUNTER} {n} {lens}"
            if phase == "timing":
                shell = f"/usr/bin/time -v sh -c '{shell} > /dev/null'"
            proc = subprocess.run(shell, shell=True, input=seq + "\n", capture_output=True,
                                  text=True, check=True)
            record = ({"vals": parse_counter(proc.stdout, n)} if phase == "accuracy"
                      else dict(zip(("wall_s", "max_rss_kb"), base.time_v(proc.stderr))))
        else:
            beam = 0 if cfg == "exact" else int(cfg[len("linracc_b"):])
            cmd = [str(linracc), f"-seqfile={fa}", f"-outfile={out}", f"-access_len={lens}",
                   f"-beam={beam}", "-energy=turner2004", "-c_multi=auto", "-c_hairpin=auto",
                   "-byloop"] + (["-no-fast-logsumexp"] if cfg == "exact" else [])
            if phase == "timing":
                proc = subprocess.run(["/usr/bin/time", "-v", *cmd], capture_output=True,
                                      text=True, check=True)
                record = dict(zip(("wall_s", "max_rss_kb"), base.time_v(proc.stderr)))
            else:
                subprocess.run(cmd, capture_output=True, text=True, check=True)
                record = {"vals": parse_linracc_byloop(out, n)}
    record.update(name=name, length=n, config=cfg)
    cache.write_text(json.dumps(record))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=("bprna52", "long"), required=True)
    ap.add_argument("--phase", choices=("accuracy", "timing"), required=True)
    ap.add_argument("--work-dir", type=Path, required=True)
    ap.add_argument("--jobs", type=int, default=16)
    args = ap.parse_args()
    work = args.work_dir.resolve()
    linracc = work / "bin/LinRacc"
    if not linracc.is_file():
        raise SystemExit(f"missing frozen binary {linracc}")
    (work / args.phase).mkdir(parents=True, exist_ok=True)
    (work / "tmp").mkdir(exist_ok=True)
    if args.set == "long":
        fasta = base.read_fasta_dict(base.LONG_FASTA)
        chosen = base.select_long(fasta)
    else:
        fasta = base.read_fasta_dict(base.FASTA)
        chosen = [n for n, _ in base.select_sequences(13)]
    tasks = [(n, c) for n in chosen for c in configs(args.phase)]
    if args.phase == "timing":
        for n, c in tasks:
            run_one(work, linracc, "timing", n, fasta[n], c)
        return 0
    tasks.sort(key=lambda t: -(len(fasta[t[0]]) ** (3 if t[1] == "exact" else 1)))
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(run_one, work, linracc, "accuracy", n, fasta[n], c) for n, c in tasks]
        for i, f in enumerate(as_completed(futs), 1):
            f.result()
            if i % 20 == 0:
                print(f"{i}/{len(tasks)}", file=sys.stderr, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
