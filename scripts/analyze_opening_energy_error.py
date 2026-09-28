#!/usr/bin/env python3
"""Opening-energy view of the LinearRaccess vs LinearSampling comparison.

For each configuration cached by run_linearsampling_comparison_capb.py
(52-sequence set and --set long), convert accessibility P to opening energy
E = -RT ln P (RT = 0.61633 kcal/mol, Raccess's output convention) and report,
against the exact reference:
  - |E - E_exact| in kcal/mol by window length (median, p95), over windows
    with finite exact energy; LinearSampling windows with zero sampled
    frequency have E = +inf and are counted separately (frac_infinite);
  - ranking agreement: Spearman correlation between approximate and exact
    energies across the windows of each (sequence, w), summarized as the
    median over (sequence, w) pairs.
Only sequences longer than the smallest beam compared are informative
(for n <= b LinearRaccess is exact); pass --min-length accordingly.
"""
from __future__ import annotations

import argparse
import glob
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

RT = 0.61633008
REPO_ROOT = Path(__file__).resolve().parent.parent.parent


def energy(p: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore"):
        return -RT * np.log(np.clip(p, 0.0, None))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work-dir", type=Path, required=True)
    ap.add_argument("--min-length", type=int, default=0)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    runs: dict[str, dict[str, dict]] = defaultdict(dict)
    for f in glob.glob(str(args.work_dir / "accuracy" / "*.json")):
        d = json.load(open(f))
        runs[d["name"]][d["config"]] = d
    names = [n for n, r in runs.items() if r["exact"]["length"] > args.min_length]
    configs = sorted(c for c in runs[names[0]] if c != "exact")

    summary = {"n_sequences": len(names), "min_length": args.min_length, "configs": {}}
    for c in configs:
        per_w = defaultdict(lambda: {"err": [], "inf": 0, "n": 0, "rho": []})
        for n in names:
            ex, ap_ = runs[n]["exact"]["vals"], runs[n][c]["vals"]
            for w in ex:
                e_ex, e_ap = energy(np.array(ex[w])), energy(np.array(ap_[w]))
                ok = np.isfinite(e_ex)
                acc = per_w[w]
                acc["n"] += int(ok.sum())
                acc["inf"] += int((ok & ~np.isfinite(e_ap)).sum())
                fin = ok & np.isfinite(e_ap)
                acc["err"].extend(np.abs(e_ap[fin] - e_ex[fin]).tolist())
                if ok.sum() > 2:
                    # +inf ranks last, as the least accessible windows.
                    a = np.where(np.isfinite(e_ap), e_ap, np.inf)[ok]
                    rho = spearmanr(a, e_ex[ok]).correlation
                    if not math.isnan(rho):
                        acc["rho"].append(rho)
        out = {}
        for w, acc in sorted(per_w.items(), key=lambda kv: int(kv[0])):
            err = np.array(acc["err"])
            out[w] = {"median_abs_kcal": float(np.median(err)),
                      "p95_abs_kcal": float(np.percentile(err, 95)),
                      "frac_infinite": acc["inf"] / acc["n"],
                      "median_spearman": float(np.median(acc["rho"])),
                      "min_spearman": float(np.min(acc["rho"]))}
        summary["configs"][c] = out
    text = json.dumps(summary, indent=1)
    if args.out:
        args.out.write_text(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
