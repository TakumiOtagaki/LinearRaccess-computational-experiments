#!/usr/bin/env python3
"""Does LinearRaccess's beam=100 pruning error depend on the energy model?

Table 1 / S3.2 (Raccess-compatible Turner 1999, -engine=raccess) reported a
median relative error of 0.6% on low-probability cells, but under the
default Turner 2004 (-engine=lincapr) the same measure is ~42%
(run_linracc_error_decomposition.py). Same binary (build-raccess, which
supports both engines), same 52 bpRNA sequences, for each energy model:
beam=0/caps=len (reference) vs beam=100/caps=30.
"""
from __future__ import annotations

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
WORK_DIR = REPO_ROOT / "benchmark/results/.work/linracc_energy_model_beam_check"
ACCESS_LENS = [1, 3, 7, 10]
MODELS = {
    "raccessT99": ["-engine=raccess", "-energy=turner1999"],
    "lincaprT99": ["-engine=lincapr", "-energy=turner1999"],
    "lincaprT04": ["-engine=lincapr", "-energy=turner2004"],
}
BEAMS = {"b0_caplen": (0, "len"), "b100_cap30": (100, "30")}


def run(name, seq, model, btag):
    cache = WORK_DIR / f"{name}__{model}__{btag}.json"
    if cache.exists():
        return
    beam, cap = BEAMS[btag]
    capv = str(len(seq)) if cap == "len" else cap
    with tempfile.TemporaryDirectory() as td:
        fa, out = Path(td) / "s.fa", Path(td) / "s.out"
        fa.write_text(f">{name}\n{seq}\n")
        cmd = [str(LINRACC), f"-seqfile={fa}", f"-outfile={out}",
               f"-access_len={','.join(map(str, ACCESS_LENS))}", f"-beam={beam}",
               f"-c_hairpin={capv}", f"-c_multi={capv}", "-probabilities", *MODELS[model]]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        parsed = parse_linracc(out)[name]
    cache.write_text(json.dumps(
        {str(w): [parsed[p][w] for p in range(len(seq) - w + 1)] for w in ACCESS_LENS}))


def main():
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    fasta = read_fasta_dict(FASTA)
    chosen = select_sequences(13)
    tasks = [(n, m, b) for n, _ in chosen for m in MODELS for b in BEAMS]
    tasks.sort(key=lambda t: -len(fasta[t[0]]))
    with ThreadPoolExecutor(max_workers=16) as ex:
        futs = [ex.submit(run, n, fasta[n], m, b) for n, m, b in tasks]
        for i, f in enumerate(as_completed(futs), 1):
            f.result()
            if i % 20 == 0:
                print(f"{i}/{len(tasks)}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
