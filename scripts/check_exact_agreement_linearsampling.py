#!/usr/bin/env python3
"""Check that LinearRaccess's exact mode (beam=0, exact log-sum-exp, Turner
2004) agrees with LinearSampling's exact sampling (-b 0) within sampling
noise, justifying the former as the reference of Supplementary S5 and S7.
For every window with exact accessibility P >= 1e-3 (w in {1,3,7,10,35}),
the z-score (P_sample - P_exact) / sqrt(P(1-P)/k) is computed, both for the
total and for each loop context.  Writes results/exact_agreement_linearsampling.json.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "benchmark/scripts"))
import run_byloop_comparison as bl  # noqa: E402
from run_linearsampling_accuracy_vs_k import read_fasta_dict  # noqa: E402

LINRACC = REPO / "vendor/LinearRaccess/build-raccess/LinRacc"
IDS = ["bpRNA_CRW_1328", "bpRNA_CRW_1563", "bpRNA_CRW_1041"]  # 288, 476, 654 nt
K = 200_000


def main() -> int:
    fasta = read_fasta_dict(REPO / "benchmark/data/bprna_100_1000.fa")
    lens = ",".join(map(str, bl.ACCESS_LENS))
    out = {"k": K, "sequences": {}, "linracc_build": subprocess.run(
        [str(LINRACC), "--build-info"], capture_output=True, text=True).stdout.strip()}
    for sid in IDS:
        seq = fasta[sid]
        with tempfile.TemporaryDirectory() as td:
            fa, o = Path(td) / "s.fa", Path(td) / "s.out"
            fa.write_text(f">{sid}\n{seq}\n")
            subprocess.run([str(LINRACC), f"-seqfile={fa}", f"-outfile={o}", f"-access_len={lens}",
                            "-beam=0", "-energy=turner2004", "-no-fast-logsumexp", "-byloop"],
                           check=True, capture_output=True)
            exact = bl.parse_linracc_byloop(o, len(seq))
        proc = subprocess.run(f"{bl.LINSAMP} -k {K} -b 0 | tail -n +2 | {bl.COUNTER} {len(seq)} {lens}",
                              shell=True, input=seq + "\n", capture_output=True, text=True, check=True)
        samp = bl.parse_counter(proc.stdout, len(seq))
        rec = {}
        for t in list(bl.TYPES) + ["total"]:
            zs = []
            for w in map(str, bl.ACCESS_LENS):
                for i in range(len(seq) - int(w) + 1):
                    pe = sum(exact[c][w][i] for c in bl.TYPES) if t == "total" else exact[t][w][i]
                    ps = sum(samp[c][w][i] for c in bl.TYPES) if t == "total" else samp[t][w][i]
                    if pe >= 1e-3 and pe < 1:
                        zs.append((ps - pe) / math.sqrt(pe * (1 - pe) / K))
            rec[t] = {"n_windows": len(zs), "max_abs_z": max(map(abs, zs)) if zs else None,
                      "n_abs_z_over_4": sum(abs(z) > 4 for z in zs)}
        out["sequences"][sid] = {"length": len(seq), **rec}
    path = REPO / "benchmark/results/exact_agreement_linearsampling.json"
    path.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
