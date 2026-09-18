#!/usr/bin/env python3
"""Regenerate Supplementary Figure S1 (per-window error distribution).

Recreated because the original figure had no retained source script and used
a color for Raccess Window (green) that collides with Figure 2/3's
convention (green = Raccess Exact). This script matches Figure 2/3's exact
color codes and font settings so all figures encode method identity
consistently across the main text and Supplementary Information.

Usage:
    python3 benchmark/scripts/plot_supp_s1_error_distribution.py

Reads benchmark/results/bprna_accuracy_comparison/{exact,linear_b100,window_w100}
for the 12,872 sequences in benchmark/data/bprna_100_1000_list.txt, and
writes manuscript/figures/figureS1_error_distribution.pdf.
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

COLORS = {
    "LinearRaccess": "#d62728",              # Red -- matches Figure 2/3
    "Raccess Window (W=100)": "#1f77b4",      # Blue -- matches Figure 2/3
}

BINS = [(100, 300), (300, 500), (500, 700), (700, 1001)]  # last bin closed at 1000
BIN_LABELS = ["100-300", "300-500", "500-700", "700-1000"]

EXACT_DIR = REPO_ROOT / "benchmark/results/bprna_accuracy_comparison/exact"
B100_DIR = REPO_ROOT / "benchmark/results/bprna_accuracy_comparison/linear_b100"
WINDOW_DIR = REPO_ROOT / "benchmark/results/bprna_accuracy_comparison/window_w100"
SEQ_LIST = REPO_ROOT / "benchmark/data/bprna_100_1000_list.txt"
OUTPUT = REPO_ROOT / "manuscript/figures/figureS1_error_distribution.pdf"


def load_list() -> list[tuple[str, int]]:
    seqs = []
    for line in SEQ_LIST.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        seq_id, length = line.split()
        seqs.append((seq_id, int(length)))
    return seqs


def bin_for_length(length: int) -> int | None:
    for i, (lo, hi) in enumerate(BINS):
        if lo <= length < hi:
            return i
    return None


def all_errors(exact_map: dict, approx_map: dict) -> list[float]:
    errors = []
    for pos, lens in exact_map.items():
        arow = approx_map.get(pos)
        if not arow:
            continue
        for length, exact_energy in lens.items():
            if length not in arow:
                continue
            p_exact = probability_from_energy(exact_energy)
            errors.append(abs(p_exact - arow[length]))
    return errors


def main() -> int:
    seqs = load_list()
    errors_by_bin = {method: [[] for _ in BINS] for method in ("LinearRaccess", "Raccess Window (W=100)")}

    for seq_id, length in seqs:
        b = bin_for_length(length)
        if b is None:
            continue
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
        errors_by_bin["LinearRaccess"][b].extend(all_errors(exact_map, lr_map))
        errors_by_bin["Raccess Window (W=100)"][b].extend(all_errors(exact_map, window_prob_map))

    fig, ax = plt.subplots(figsize=(10, 6))
    width = 0.35
    positions = np.arange(len(BINS))

    # Exact-zero errors are real, valid data points (Window matches Exact
    # exactly when no base pair spans >100nt in that window), not missing
    # data -- they must stay in the median calculation to match Table 1's
    # own (linear-space) medians. A log-scale violin plot cannot show an
    # exact zero, so we floor it far below the smallest real nonzero error
    # observed (~1e-27) rather than dropping it, which would otherwise bias
    # the displayed median upward by excluding the lowest values.
    ZERO_FLOOR = 1e-30

    for method_idx, method in enumerate(("LinearRaccess", "Raccess Window (W=100)")):
        offset = (method_idx - 0.5) * width
        data = []
        for b in range(len(BINS)):
            vals = np.array(errors_by_bin[method][b])
            vals = np.maximum(vals, ZERO_FLOOR)
            log_vals = np.log10(vals)
            data.append(log_vals)
        parts = ax.violinplot(
            data,
            positions=positions + offset,
            widths=width * 0.9,
            showmeans=False,
            showmedians=True,
            showextrema=True,
        )
        for pc in parts["bodies"]:
            pc.set_facecolor(COLORS[method])
            pc.set_alpha(0.6)
            pc.set_edgecolor(COLORS[method])
        for key in ("cmedians", "cmins", "cmaxes", "cbars"):
            parts[key].set_edgecolor(COLORS[method])

    ax.set_xticks(positions)
    ax.set_xticklabels(BIN_LABELS)
    ax.set_xlabel("Sequence length range (nt)", fontsize=16)
    ax.set_ylabel(r"$\log_{10}$(absolute error)", fontsize=16)
    ax.set_title(
        "Per-window absolute error distribution\nLinearRaccess (beam=100) vs Raccess Window (W=100)",
        fontsize=15,
    )
    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=COLORS["LinearRaccess"], alpha=0.6, label="LinearRaccess (beam=100)"),
        plt.Rectangle((0, 0), 1, 1, facecolor=COLORS["Raccess Window (W=100)"], alpha=0.6, label="Raccess Window (W=100)"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=12)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUTPUT)
    print(f"Wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
