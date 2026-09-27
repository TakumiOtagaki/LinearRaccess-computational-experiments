#!/usr/bin/env python3
"""LinearRaccess vs LinearSampling (docs/pi_feedback_20260925.md), rerun for
the single-pass accessibility build with loop caps tied to the beam and
w in {1,3,7,10,35}.

Same 52 bpRNA sequences (13 per Table 1 length bin) and Turner 2004 energy
model as run_beam_sweep_vs_exact.py / run_exact_reference_si.py.

--phase accuracy (parallel; outputs only, no timing):
  exact        LinearRaccess beam=0 (caps auto = N), exact log-sum-exp;
               the common exact reference for both tools
  linracc_b{B} LinearRaccess beam=B, caps auto (= B), default fast LSE
  linsamp_b{B}_k{K}  LinearSampling beam=B, K samples, counted per window by
               window_freq_counter (C++)
--phase timing (strictly sequential, one process at a time; /usr/bin/time -v
  around the whole command, i.e. LinearSampling includes the sample
  counting, since that is what produces accessibility):
  the same LinearRaccess / LinearSampling configurations.

Results are cached per (sequence, configuration) under WORK_DIR so either
phase can be resumed.  Analysis is done separately from the caches.

--set long: the same comparison on sequences longer than the beam, where
pruning is active (for n <= b LinearRaccess is exact): 5 bpRNA-1m(90)
sequences per bin of 1,000-1,500 / 1,500-2,000 / 2,000-3,000 /
3,000-4,381 nt, in file order, with beams 100/200/500 for both tools
(LinearSampling at k=1e6), under benchmark/results/linearsampling_comparison_long.
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
from run_raccess_compare import parse_linracc  # noqa: E402
from run_linearsampling_accuracy_vs_k import read_fasta_dict, select_sequences  # noqa: E402

WORK_DIR = REPO_ROOT / "benchmark/results/linearsampling_comparison_capb"
LINRACC = WORK_DIR / "bin/LinRacc"
LINSAMP = REPO_ROOT / "vendor/LinearSampling/linearsampling"
COUNTER = REPO_ROOT / "benchmark/scripts/window_freq_counter"
FASTA = REPO_ROOT / "benchmark/data/bprna_100_1000.fa"
ACCESS_LENS = [1, 3, 7, 10, 35]
LINRACC_BEAMS = [100, 200, 300, 500]
LINSAMP_ACCURACY = [(b, 1_000_000) for b in (100, 200, 300, 500)] + [(100, 50_000)]
LINSAMP_TIMING = [(100, 50_000), (100, 1_000_000), (500, 50_000), (500, 1_000_000)]

LONG_FASTA = REPO_ROOT / "benchmark/data/bpRNA_1m_90_filtered_clean.fa"
LONG_BINS = [(1000, 1500), (1500, 2000), (2000, 3000), (3000, 4382)]
LONG_PER_BIN = 5


def use_long_set() -> None:
    global WORK_DIR, LINRACC, FASTA, LINRACC_BEAMS, LINSAMP_ACCURACY, LINSAMP_TIMING
    WORK_DIR = REPO_ROOT / "benchmark/results/linearsampling_comparison_long"
    LINRACC = WORK_DIR / "bin/LinRacc"
    FASTA = LONG_FASTA
    LINRACC_BEAMS = [100, 200, 500]
    LINSAMP_ACCURACY = [(b, 1_000_000) for b in (100, 200, 500)]
    LINSAMP_TIMING = list(LINSAMP_ACCURACY)


def select_long(fasta: dict[str, str]) -> list[str]:
    chosen = [[] for _ in LONG_BINS]
    for name, seq in fasta.items():
        for i, (lo, hi) in enumerate(LONG_BINS):
            if lo <= len(seq) < hi and len(chosen[i]) < LONG_PER_BIN:
                chosen[i].append(name)
    return [n for group in chosen for n in group]


def linracc_cmd(fa: Path, out: Path, beam: int, exact_lse: bool) -> list[str]:
    cmd = [str(LINRACC), f"-seqfile={fa}", f"-outfile={out}",
           f"-access_len={','.join(map(str, ACCESS_LENS))}", f"-beam={beam}",
           "-energy=turner2004", "-c_multi=auto", "-c_hairpin=auto", "-probabilities"]
    return cmd + (["-no-fast-logsumexp"] if exact_lse else [])


def linsamp_shell(length: int, beam: int, k: int) -> str:
    return (f"{LINSAMP} -k {k} -b {beam} | tail -n +2 | "
            f"{COUNTER} {length} {','.join(map(str, ACCESS_LENS))}")


def parse_counter(stdout: str, length: int) -> dict[str, list[float]]:
    res = {str(w): [0.0] * (length - w + 1) for w in ACCESS_LENS}
    for line in stdout.splitlines():
        w, pos, f = line.split()
        res[w][int(pos)] = float(f)
    return res


def time_v(stderr: str) -> tuple[float, int]:
    m = re.search(r"Elapsed \(wall clock\) time.*: (?:(\d+):)?(\d+):([\d.]+)", stderr)
    wall = int(m.group(1) or 0) * 3600 + int(m.group(2)) * 60 + float(m.group(3))
    rss = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", stderr).group(1))
    return wall, rss


def configs(phase: str) -> list[str]:
    if phase == "accuracy":
        return (["exact"] + [f"linracc_b{b}" for b in LINRACC_BEAMS]
                + [f"linsamp_b{b}_k{k}" for b, k in LINSAMP_ACCURACY])
    return [f"linracc_b{b}" for b in LINRACC_BEAMS] + [f"linsamp_b{b}_k{k}" for b, k in LINSAMP_TIMING]


def run_one(phase: str, name: str, seq: str, cfg: str) -> None:
    cache = WORK_DIR / phase / f"{name}__{cfg}.json"
    if cache.exists():
        return
    with tempfile.TemporaryDirectory(dir=WORK_DIR / "tmp") as td:
        fa, out = Path(td) / "s.fa", Path(td) / "s.out"
        fa.write_text(f">{name}\n{seq}\n")
        if cfg.startswith("linsamp"):
            beam, k = map(int, re.fullmatch(r"linsamp_b(\d+)_k(\d+)", cfg).groups())
            shell = linsamp_shell(len(seq), beam, k)
            if phase == "timing":
                shell = f"/usr/bin/time -v sh -c '{shell} > /dev/null'"
            proc = subprocess.run(shell, shell=True, input=seq + "\n", capture_output=True,
                                  text=True, check=True)
            record = ({"vals": parse_counter(proc.stdout, len(seq))} if phase == "accuracy"
                      else dict(zip(("wall_s", "max_rss_kb"), time_v(proc.stderr))))
        else:
            beam = 0 if cfg == "exact" else int(cfg[len("linracc_b"):])
            cmd = linracc_cmd(fa, out, beam, exact_lse=(cfg == "exact"))
            if phase == "timing":
                proc = subprocess.run(["/usr/bin/time", "-v", *cmd], capture_output=True,
                                      text=True, check=True)
                record = dict(zip(("wall_s", "max_rss_kb"), time_v(proc.stderr)))
            else:
                subprocess.run(cmd, capture_output=True, text=True, check=True)
                parsed = parse_linracc(out)[name]
                record = {"vals": {str(w): [parsed[p][w] for p in range(len(seq) - w + 1)]
                                   for w in ACCESS_LENS}}
    record.update(name=name, length=len(seq), config=cfg)
    cache.write_text(json.dumps(record))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=("accuracy", "timing"), required=True)
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--set", choices=("bprna52", "long"), default="bprna52")
    args = parser.parse_args()
    if args.set == "long":
        use_long_set()
    if not LINRACC.is_file():
        raise SystemExit(f"missing frozen binary {LINRACC}")
    (WORK_DIR / args.phase).mkdir(parents=True, exist_ok=True)
    (WORK_DIR / "tmp").mkdir(exist_ok=True)
    fasta = read_fasta_dict(FASTA)
    chosen = select_long(fasta) if args.set == "long" else [name for name, _ in select_sequences(13)]
    tasks = [(n, c) for n in chosen for c in configs(args.phase)]
    if args.phase == "timing":
        for i, (n, c) in enumerate(tasks, 1):
            run_one("timing", n, fasta[n], c)
            if i % 20 == 0:
                print(f"{i}/{len(tasks)}", file=sys.stderr, flush=True)
        return 0
    # Heaviest (exact, long) first so the tail of the run is short.
    tasks.sort(key=lambda t: -(len(fasta[t[0]]) ** (3 if t[1] == "exact" else 1)))
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(run_one, "accuracy", n, fasta[n], c) for n, c in tasks]
        for i, f in enumerate(as_completed(futs), 1):
            f.result()
            if i % 20 == 0:
                print(f"{i}/{len(tasks)}", file=sys.stderr, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
