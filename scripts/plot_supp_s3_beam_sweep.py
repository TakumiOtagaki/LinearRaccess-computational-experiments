#!/usr/bin/env python3
"""Supplementary Figure S3: LinearRaccess error vs beam width on the 1,770
bpRNA sequences of 300-1,000 nt, against Raccess Exact (Turner 1999,
w in {1,3,7,10,35}, loop caps tied to the beam).

Reads the aggregates written by accuracy_aggregate for LinearRaccess 0155fa4:
benchmark/results/bprna_accuracy_dd_summary/sweep_b{50,100,200,300,500}.json
(A) median absolute error, (B) median opening-energy error, overall and per
length bin; length bins are ordered, so they share one hue light -> dark.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
SUMMARY = REPO_ROOT / "benchmark/results/bprna_accuracy_dd_summary"
OUTPUT = REPO_ROOT / "manuscript/figures/figureS3_beam_sweep.pdf"
BEAMS = [50, 100, 200, 300, 500]
SERIES = [  # (bin key, label, color, marker)
    ("300-500", "300--500 nt", "#9ecae1", "o"),
    ("500-700", "500--700 nt", "#4292c6", "s"),
    ("700-1000", "700--1,000 nt", "#08519c", "^"),
    ("pooled", "all (1,770)", "#222222", "D"),
]
FLOOR = 1e-10

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"]
plt.rcParams["font.size"] = 12


def main() -> int:
    data = {b: json.loads((SUMMARY / f"sweep_b{b}.json").read_text())["bins"] for b in BEAMS}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    for ax, key, ylabel, panel in [
        (axes[0], "median_abs", "Median absolute error", "A"),
        (axes[1], "median_dG_kcal", "Median opening-energy error (kcal/mol)", "B"),
    ]:
        for bin_key, label, color, marker in SERIES:
            ys = [max(data[b][bin_key]["all"][key], FLOOR) for b in BEAMS]
            ax.plot(BEAMS, ys, marker=marker, color=color, linewidth=2, markersize=8,
                    label=label.replace("--", "–"))
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks(BEAMS)
        ax.set_xticklabels([str(b) for b in BEAMS])
        ax.set_xlabel("Beam width $b$")
        ax.set_ylabel(ylabel)
        ax.grid(True, which="major", color="#dddddd", linewidth=0.6)
        ax.text(-0.13, 1.04, f"({panel})", transform=ax.transAxes, fontsize=14, fontweight="bold")
    axes[0].legend(frameon=False, fontsize=10, loc="lower left")
    fig.tight_layout()
    fig.savefig(OUTPUT, bbox_inches="tight")
    fig.savefig(OUTPUT.with_suffix(".png"), dpi=300, bbox_inches="tight")
    print(f"Wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
