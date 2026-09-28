#!/usr/bin/env python3
"""Regenerate the Results \S3.2 scatter-plot figure: LinearRaccess and Raccess
Window accessibility values plotted against Raccess Exact (ground truth), one
panel per method.

Reads benchmark/results/bprna_accuracy_capb/{exact,window_w100} and
bprna_accuracy_dd/linear_b200 (LinearRaccess 0155fa4, beam=200, w in {1,3,7,10,35})
for the 12,872 sequences in benchmark/data/bprna_100_1000_list.txt (the same
evidence bundle and sequence set as Table 1 and Figure S1), and writes
manuscript/figures/figure_scatter_accuracy.pdf.

Usage:
    python3 benchmark/scripts/plot_scatter_accuracy.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(SCRIPT_DIR))

from run_raccess_compare import parse_raccess, parse_linracc, probability_from_energy

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"]
plt.rcParams["font.size"] = 14

EXACT_DIR = REPO_ROOT / "benchmark/results/bprna_accuracy_capb/exact"
B100_DIR = REPO_ROOT / "benchmark/results/bprna_accuracy_dd/linear_b200"
WINDOW_DIR = REPO_ROOT / "benchmark/results/bprna_accuracy_capb/window_w100"
SEQ_LIST = REPO_ROOT / "benchmark/data/bprna_100_1000_list.txt"
OUTPUT = REPO_ROOT / "manuscript/figures/figure_scatter_accuracy.pdf"

# Floor for the log-scale axes; matches the numerical-tolerance floor used
# elsewhere in this paper's evidence bundles.
FLOOR = 1e-12


def load_list() -> list[tuple[str, int]]:
    seqs = []
    for line in SEQ_LIST.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        seq_id, length = line.split()
        seqs.append((seq_id, int(length)))
    return seqs


def paired_values(exact_map: dict, approx_map: dict) -> tuple[list[float], list[float]]:
    exact_vals: list[float] = []
    approx_vals: list[float] = []
    for pos, lens in exact_map.items():
        arow = approx_map.get(pos)
        if not arow:
            continue
        for length, exact_energy in lens.items():
            if length not in arow:
                continue
            exact_vals.append(probability_from_energy(exact_energy))
            approx_vals.append(arow[length])
    return exact_vals, approx_vals


def main() -> int:
    seqs = load_list()

    exact_lr: list[float] = []
    lr_vals: list[float] = []
    exact_win: list[float] = []
    win_vals: list[float] = []

    for seq_id, _length in seqs:
        fname = seq_id + ".txt"
        exact_path = EXACT_DIR / fname
        b100_path = B100_DIR / fname
        window_path = WINDOW_DIR / fname
        if not (exact_path.exists() and b100_path.exists() and window_path.exists()):
            continue
        exact_map = parse_raccess(exact_path).get(seq_id, {})
        lr_map = parse_linracc(b100_path).get(seq_id, {})
        window_map = parse_raccess(window_path).get(seq_id, {})
        window_prob_map = {
            pos: {length: probability_from_energy(e) for length, e in lens.items()}
            for pos, lens in window_map.items()
        }

        e1, a1 = paired_values(exact_map, lr_map)
        exact_lr.extend(e1)
        lr_vals.extend(a1)

        e2, a2 = paired_values(exact_map, window_prob_map)
        exact_win.extend(e2)
        win_vals.extend(a2)

    print(f"LinearRaccess vs Exact: {len(exact_lr):,} cells")
    print(f"Raccess Window vs Exact: {len(exact_win):,} cells")

    exact_lr_arr = np.maximum(np.array(exact_lr), FLOOR)
    lr_arr = np.maximum(np.array(lr_vals), FLOOR)
    exact_win_arr = np.maximum(np.array(exact_win), FLOOR)
    win_arr = np.maximum(np.array(win_vals), FLOOR)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    panels = [
        (axes[0], exact_lr_arr, lr_arr, "LinearRaccess (beam=200)", "A"),
        (axes[1], exact_win_arr, win_arr, "Raccess Window (W=100)", "B"),
    ]

    for ax, exact_arr, approx_arr, label, tag in panels:
        hb = ax.hexbin(
            exact_arr,
            approx_arr,
            gridsize=60,
            bins="log",
            xscale="log",
            yscale="log",
            cmap="viridis",
            mincnt=1,
        )
        lims = [FLOOR, 1.0]
        ax.plot(lims, lims, color="black", linewidth=1.2, linestyle="--", alpha=0.8, label="$y=x$")
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel("Raccess Exact accessibility", fontsize=16)
        ax.set_ylabel(f"{label} accessibility", fontsize=16)
        ax.set_title(f"({tag})", fontsize=18, fontweight="bold", loc="left")
        ax.legend(loc="upper left", fontsize=12)
        cb = fig.colorbar(hb, ax=ax)
        cb.set_label("Count (log scale)", fontsize=12)

    fig.tight_layout()
    fig.savefig(OUTPUT, bbox_inches="tight")
    fig.savefig(OUTPUT.with_suffix(".png"), dpi=300, bbox_inches="tight")
    print(f"Figure saved to: {OUTPUT}")
    print(f"Figure saved to: {OUTPUT.with_suffix('.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
