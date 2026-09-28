#!/usr/bin/env python3
"""Supplementary S7 table from run_byloop_comparison.py caches: for each loop
type (E/H/B/I/M) and configuration, over windows with exact per-type
accessibility P^c >= 1e-8: median relative error, and the percentage whose
estimate is exactly zero; plus the share of all accessibility mass per type.

  analyze_byloop_comparison.py --work-dir DIR [--min-length N] --out JSON
"""
from __future__ import annotations

import argparse
import glob
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

TYPES = "EHBIM"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work-dir", type=Path, required=True)
    ap.add_argument("--min-length", type=int, default=0)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    runs = defaultdict(dict)
    for f in glob.glob(str(args.work_dir / "accuracy/*.json")):
        d = json.load(open(f))
        runs[d["name"]][d["config"]] = d
    names = sorted(n for n in runs if runs[n]["exact"]["length"] > args.min_length)
    out = {"n_sequences": len(names), "min_length": args.min_length, "configs": {}, "mass_share": {}}
    mass = {t: sum(float(np.sum(v)) for n in names for v in runs[n]["exact"]["vals"][t].values()) for t in TYPES}
    tot = sum(mass.values())
    out["mass_share"] = {t: mass[t] / tot for t in TYPES}
    for c in sorted(runs[names[0]]):
        if c == "exact":
            continue
        rec = {}
        for t in TYPES:
            rel, zero, n = [], 0, 0
            for nm in names:
                for w, ev in runs[nm]["exact"]["vals"][t].items():
                    e = np.array(ev); a = np.array(runs[nm][c]["vals"][t][w])
                    ok = e >= 1e-8
                    rel.append(np.abs(a[ok] - e[ok]) / e[ok]); zero += int(np.sum(a[ok] <= 0)); n += int(ok.sum())
            rel = np.concatenate(rel)
            rec[t] = {"n_windows": n, "median_rel": float(np.median(rel)), "pct_zero": 100.0 * zero / n}
        out["configs"][c] = rec
    args.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
