#!/usr/bin/env python3
"""Regenerate Supplementary Figure S2 (effect of beam width on hard sequences).

Plots the distribution of PER-SEQUENCE MAXIMUM absolute error (one point per
affected sequence), matching the main-text/Supplementary S1 prose ("median
per-sequence maximum error drops from 0.785 (beam=100) to 0.032
(beam=200)"). An earlier version of this script pooled all per-window errors
instead, which is a different statistic (median error ~1e-3, not directly
comparable to the per-sequence maxima quoted in prose) and did not match the
number readers would try to read off the plot; fixed after review caught the
mismatch (see docs/review_triage_20260918.md).

Usage:
    python3 benchmark/scripts/plot_supp_s2_beam_width_effect.py

Reads benchmark/results/bprna_accuracy_capb/exact and
bprna_accuracy_dd/{linear_b100,linear_b200} (LinearRaccess 0155fa4, w in {1,3,7,10,35})
for the 791 sequences with max absolute error > 0.5 at beam=100 (identified
from the full 12,872-sequence set), and writes
manuscript/figures/figureS2_beam_width_effect.pdf.
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
    "beam=100": "#d62728",  # Red -- matches Figure 2/3/S1
    "beam=200": "#1f77b4",  # Blue -- matches Figure 2/3/S1
}

EXACT_DIR = REPO_ROOT / "benchmark/results/bprna_accuracy_capb/exact"
B100_DIR = REPO_ROOT / "benchmark/results/bprna_accuracy_dd/linear_b100"
B200_DIR = REPO_ROOT / "benchmark/results/bprna_accuracy_dd/linear_b200"
SEQ_LIST = REPO_ROOT / "benchmark/data/bprna_100_1000_list.txt"
OUTPUT = REPO_ROOT / "manuscript/figures/figureS2_beam_width_effect.pdf"

MAX_ERROR_THRESHOLD = 0.5


def load_ids() -> list[str]:
    ids = []
    for line in SEQ_LIST.read_text().splitlines():
        line = line.strip()
        if line:
            ids.append(line.split()[0])
    return ids


def all_errors_and_max(exact_map: dict, approx_map: dict) -> tuple[list[float], float]:
    errors = []
    max_err = 0.0
    for pos, lens in exact_map.items():
        arow = approx_map.get(pos)
        if not arow:
            continue
        for length, exact_energy in lens.items():
            if length not in arow:
                continue
            p_exact = probability_from_energy(exact_energy)
            err = abs(p_exact - arow[length])
            errors.append(err)
            if err > max_err:
                max_err = err
    return errors, max_err


def main() -> int:
    seq_ids = load_ids()

    hard_ids = []
    for seq_id in seq_ids:
        fname = seq_id + ".txt"
        exact_path = EXACT_DIR / fname
        b100_path = B100_DIR / fname
        if not (exact_path.exists() and b100_path.exists()):
            continue
        exact_map = parse_raccess(exact_path).get(seq_id, {})
        b100_map = parse_linracc(b100_path).get(seq_id, {})
        _, max_err = all_errors_and_max(exact_map, b100_map)
        if max_err > MAX_ERROR_THRESHOLD:
            hard_ids.append(seq_id)

    print(f"n hard sequences (max error > {MAX_ERROR_THRESHOLD} at beam=100): {len(hard_ids)}")

    max_errors_beam100: list[float] = []
    max_errors_beam200: list[float] = []
    for seq_id in hard_ids:
        fname = seq_id + ".txt"
        exact_map = parse_raccess(EXACT_DIR / fname).get(seq_id, {})
        b100_map = parse_linracc(B100_DIR / fname).get(seq_id, {})
        _, max_err100 = all_errors_and_max(exact_map, b100_map)
        max_errors_beam100.append(max_err100)

        b200_path = B200_DIR / fname
        if b200_path.exists():
            b200_map = parse_linracc(b200_path).get(seq_id, {})
            _, max_err200 = all_errors_and_max(exact_map, b200_map)
            max_errors_beam200.append(max_err200)

    print(f"median per-sequence max error at beam=100: {np.median(max_errors_beam100):.3f}")
    print(f"median per-sequence max error at beam=200: {np.median(max_errors_beam200):.3f}")

    fig, ax = plt.subplots(figsize=(6, 6))
    labels = ["beam=100", "beam=200"]
    datasets = [np.array(max_errors_beam100), np.array(max_errors_beam200)]
    positions = [1, 2]

    # Floor exact-zero errors instead of dropping them (see
    # plot_supp_s1_error_distribution.py for why: they are real data points
    # and dropping them biases the displayed median upward).
    ZERO_FLOOR = 1e-30

    for pos, label, data in zip(positions, labels, datasets):
        data = np.maximum(data, ZERO_FLOOR)
        log_data = np.log10(data)
        parts = ax.violinplot([log_data], positions=[pos], widths=0.8,
                               showmeans=False, showmedians=True, showextrema=True)
        for pc in parts["bodies"]:
            pc.set_facecolor(COLORS[label])
            pc.set_alpha(0.6)
            pc.set_edgecolor(COLORS[label])
        for key in ("cmedians", "cmins", "cmaxes", "cbars"):
            parts[key].set_edgecolor(COLORS[label])

    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    ax.set_ylabel(r"$\log_{10}$(per-sequence maximum absolute error)", fontsize=15)
    ax.set_title(
        f"Effect of beam width on the {len(hard_ids)} sequences\nwith max error > {MAX_ERROR_THRESHOLD} at beam=100",
        fontsize=15,
    )
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUTPUT)
    print(f"Wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
