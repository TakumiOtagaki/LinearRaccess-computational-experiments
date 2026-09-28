#!/usr/bin/env python3
"""Summary tables for the LinearRaccess vs LinearSampling comparison
(Supplementary S5) from run_linearsampling_comparison_capb.py caches:
median relative error by exact-probability decade (all windows,
w in {1,3,7,10,35}), overall median/p95 relative error over windows with
P >= 1e-8, sequences with some |dP| > 0.5, and single-core time / peak RSS.

  analyze_linearsampling_comparison.py --work-dir DIR [--min-length N] --out JSON
"""
from __future__ import annotations

import argparse
import glob
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

DECADES = [(1e-1, 1.01), (1e-2, 1e-1), (1e-3, 1e-2), (1e-4, 1e-3), (1e-6, 1e-4), (1e-8, 1e-6)]


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
    timing = defaultdict(lambda: {"wall_s": 0.0, "max_rss_kb": 0, "n": 0})
    for f in glob.glob(str(args.work_dir / "timing/*.json")):
        d = json.load(open(f))
        if d["length"] > args.min_length:
            t = timing[d["config"]]
            t["wall_s"] += d["wall_s"]; t["max_rss_kb"] = max(t["max_rss_kb"], d["max_rss_kb"]); t["n"] += 1
    out = {"n_sequences": len(names), "min_length": args.min_length,
           "lengths": [runs[n]["exact"]["length"] for n in names],
           "decades": [f"[{lo:g},{hi:g})" for lo, hi in DECADES], "configs": {}}
    for c in sorted(runs[names[0]]):
        if c == "exact":
            continue
        by_dec = [[] for _ in DECADES]; allr = []; fail = 0
        for n in names:
            ex, ap_ = runs[n]["exact"]["vals"], runs[n][c]["vals"]
            worst = 0.0
            for w in ex:
                e, a = np.array(ex[w]), np.array(ap_[w])
                worst = max(worst, float(np.max(np.abs(a - e))))
                ok = e >= 1e-8
                r = np.abs(a[ok] - e[ok]) / e[ok]
                allr.append(r)
                for i, (lo, hi) in enumerate(DECADES):
                    sel = (e[ok] >= lo) & (e[ok] < hi)
                    by_dec[i].append(r[sel])
            fail += worst > 0.5
        allr = np.concatenate(allr)
        rec = {"median_rel_by_decade": [float(np.median(np.concatenate(x))) for x in by_dec],
               "median_rel_all": float(np.median(allr)), "p95_rel_all": float(np.percentile(allr, 95)),
               "n_seq_fail": int(fail)}
        if c in timing:
            rec.update(total_wall_s=timing[c]["wall_s"], max_rss_mb=timing[c]["max_rss_kb"] / 1024,
                       n_timed=timing[c]["n"])
        out["configs"][c] = rec
    args.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
