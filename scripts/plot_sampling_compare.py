#!/usr/bin/env python3
"""Regenerate the sampling-variance comparison figure from frozen bundles.

Reads one or more run_sampling_compare.py result bundles and plots:
  A) observed vs predicted relative RMS deviation of sampling estimates,
     binned by reference accessibility p_ref (points = (S, w, p-bin));
  B) observed zero-estimate fraction vs S*p_ref against the Bernoulli
     prediction (1-p)^S ~ exp(-S*p).
Used for the composite-figure panel D of the paper; every number comes from
the reproducible bundle TSVs.
"""

from __future__ import annotations

import argparse
import ast
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BINS = [1e-3, 1e-2, 1e-1, 0.5, 1.0]
BIN_MID = [0.5 * (b if i == 0 else BINS[i - 1]) + 0.5 * b if i else b / 100.0
           for i, b in enumerate(BINS)]
BIN_MID = [math.sqrt(a * b) for a, b in
           zip([0.0] + BINS[:-1], BINS)]


def load_aggregate(bundle: Path) -> list[dict]:
    with (bundle / "aggregate.tsv").open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_raw(bundle: Path) -> list[dict]:
    with (bundle / "raw.tsv").open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", "-b", action="append", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    aggregates: list[dict] = []
    raw_rows: list[dict] = []
    for bundle in args.bundle:
        if not (bundle / "aggregate.tsv").is_file():
            print(f"bundle missing aggregate.tsv: {bundle}")
            return 1
        aggregates.extend(load_aggregate(bundle))
        if (bundle / "raw.tsv").is_file():
            raw_rows.extend(load_raw(bundle))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    ax = axes[0]
    series = {}
    for row in aggregates:
        S = int(row["S"])
        pred = float(row["predicted_rel_sd"])
        obs = float(row["rmse_rel_dev"])
        key = f"S=$10^{{{math.log10(S):.0f}}}$"
        series.setdefault(key, ([], []))
        series[key][0].append(pred)
        series[key][1].append(obs)
    colors = {"S=$10^{2}$": "#4C72B0", "S=$10^{3}$": "#55A868", "S=$10^{4}$": "#C44E52",
              "S=$10^{5}$": "#8172B3"}
    markers = {"S=$10^{2}$": "o", "S=$10^{3}$": "s", "S=$10^{4}$": "^",
               "S=$10^{5}$": "D"}
    for key, (xs, ys) in sorted(series.items()):
        ax.scatter(xs, ys, s=18, alpha=0.75, label=key, color=colors.get(key),
                   marker=markers.get(key, "o"))
    ax.plot([0.005, 30.0], [0.005, 30.0], "k--", lw=0.9, label="$y=x$")
    ax.text(0.006, 1.03, "floor: all 20 estimates zero", fontsize=8,
            ha="left", va="bottom")
    ax.set_xscale("log", base=10)
    ax.set_yscale("log", base=10)
    ax.set_xlim(0.005, 30)
    ax.set_ylim(0.005, 30)
    ax.set_xlabel("predicted relative SD $\\sqrt{(1-p)/(S p)}$")
    ax.set_ylabel("observed relative RMSE (vs finite-sample reference)")
    ax.set_title("(A) Observed relative RMSE vs binomial prediction")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    name_map = {"p100": 100, "p1e3": 1000, "p1e4": 10000}
    by_s: dict[int, dict[float, list[float]]] = {}
    for row in raw_rows:
        p_ref = float(row["p_ref"])
        if p_ref <= 0:
            continue
        for key, S in name_map.items():
            estimates = row.get(key) or []
            if not isinstance(estimates, list):
                try:
                    estimates = ast.literal_eval(str(estimates))
                except (ValueError, SyntaxError):
                    estimates = []
            if not estimates:
                continue
            zf = sum(1 for e in estimates if float(e) == 0.0) / len(estimates)
            sp = S * p_ref
            if sp <= 0:
                continue
            bin_x = round(math.log10(sp) * 2) / 2
            by_s.setdefault(S, {}).setdefault(bin_x, []).append(zf)
    colors2 = {100: "#4C72B0", 1000: "#55A868", 10000: "#C44E52", 100000: "#8172B3"}
    for S in sorted(by_s):
        bins = by_s[S]
        xs = sorted(bins)
        ys = [sum(bins[x]) / len(bins[x]) for x in xs]
        ax.plot(xs, ys, marker="o", ms=4, lw=1.1,
                label=f"S=$10^{{{math.log10(S):.0f}}}$", color=colors2[S])
    sp_range = [math.log10(s) for s in (0.03, 0.1, 0.3, 1.0, 3.0, 10.0)]
    bern = [math.exp(-10.0 ** s) for s in sp_range]
    ax.plot(sp_range, bern, "k--", lw=0.9,
            label="$e^{-S p_{\\rm ref}}$ (Poisson approx.)")
    ax.set_xlabel("$\\log_{10}(S\\,p_{\\rm ref})$")
    ax.set_ylabel("mean zero-estimate fraction (bin mean)")
    ax.set_title("(B) Probability of a zero sampling estimate")
    ax.legend(frameon=False, fontsize=9, loc="upper right")

    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=300)
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())