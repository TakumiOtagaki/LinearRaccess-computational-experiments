#!/usr/bin/env python3
"""Plot Figure 3: Runtime (A) and Memory (B) as a single two-panel figure.

Item 22: the previous version placed figure2_length_vs_time_comparison.pdf
and figure3_memory.pdf side by side at 0.48\\linewidth each in the
manuscript, which shrinks their (already matplotlib-default-sized) text by
roughly 3x and makes them look cramped/inconsistent with Fig. 2's combined,
full-width panel figure (plot_scatter_accuracy.py). This script reproduces
the same two datasets in one matplotlib figure with two subplots, styled to
match plot_scatter_accuracy.py exactly (font sizes, panel-label placement,
color/marker convention), so it can be placed at \\linewidth like Fig. 2.

Usage:
    python3 benchmark/scripts/plot_figure_runtime_memory.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


def linregress_slope_intercept(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    result = stats.linregress(x, y)
    return result.slope, result.intercept

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"]
plt.rcParams["font.size"] = 14

COLORS = {
    "LinearRaccess": "#d62728",  # Red -- matches Fig. 2/S1
    "Raccess Exact": "#2ca02c",  # Green
    "Raccess Window (W=100)": "#1f77b4",  # Blue -- matches Fig. 2/S1
}
MARKERS = {
    "LinearRaccess": "o",
    "Raccess Exact": "s",
    "Raccess Window (W=100)": "^",
}

MATCHED_DIR = REPO_ROOT / "benchmark/results/bprna_matched_runtime/8714351"
RACCESS_EXACT_DIR = REPO_ROOT / "benchmark/results/raccess_exact"
OUTPUT = REPO_ROOT / "manuscript/figures/figure_runtime_memory.pdf"


def load_matched_records() -> list[dict]:
    records: list[dict] = []
    for batch_dir in sorted(MATCHED_DIR.glob("batch_*")):
        summary_file = batch_dir / "run_summary.json"
        if not summary_file.exists():
            continue
        data = json.loads(summary_file.read_text())
        records.extend(r for r in data.get("records", []) if r["linracc_ok"] and r["raccess_ok"])
    return records


def load_linearraccess_runtime() -> tuple[list[int], list[float]]:
    records = load_matched_records()
    return [r["length"] for r in records], [r["linracc_elapsed_sec"] for r in records]


def load_raccess_window_runtime() -> tuple[list[int], list[float]]:
    records = load_matched_records()
    return [r["length"] for r in records], [r["raccess_elapsed_sec"] for r in records]


def load_linearraccess_memory() -> tuple[list[int], list[float]]:
    records = [r for r in load_matched_records() if r["linracc_peak_rss_mb"] is not None]
    return [r["length"] for r in records], [r["linracc_peak_rss_mb"] for r in records]


def load_raccess_window_memory() -> tuple[list[int], list[float]]:
    records = [r for r in load_matched_records() if r["raccess_peak_rss_mb"] is not None]
    return [r["length"] for r in records], [r["raccess_peak_rss_mb"] for r in records]


def load_raccess_runtime(results_dir: Path) -> tuple[list[int], list[float]]:
    lengths: list[int] = []
    times: list[float] = []
    for json_file in sorted(results_dir.glob("*.json")):
        data = json.loads(json_file.read_text())
        if data.get("status") == "SUCCESS":
            lengths.append(data["seq_length"])
            times.append(data["elapsed_sec"])
    return lengths, times


def load_memory_data(results_dir: Path) -> tuple[list[int], list[float]]:
    lengths: list[int] = []
    memory: list[float] = []
    for summary_file in results_dir.rglob("run_summary.json"):
        data = json.loads(summary_file.read_text())
        for record in data.get("records", []):
            if record.get("status") == "SUCCESS" and record.get("peak_rss_mb") is not None:
                lengths.append(record["length"])
                memory.append(record["peak_rss_mb"])
    if not lengths:
        for json_file in results_dir.glob("*.json"):
            data = json.loads(json_file.read_text())
            if data.get("status") == "SUCCESS" and data.get("memory_kb") is not None:
                lengths.append(data["seq_length"])
                memory.append(data["memory_kb"] / 1024)
    return lengths, memory


def main() -> int:
    runtime = {
        "LinearRaccess": load_linearraccess_runtime(),
        "Raccess Exact": load_raccess_runtime(RACCESS_EXACT_DIR),
        "Raccess Window (W=100)": load_raccess_window_runtime(),
    }
    memory = {
        "LinearRaccess": load_linearraccess_memory(),
        "Raccess Exact": load_memory_data(RACCESS_EXACT_DIR),
        "Raccess Window (W=100)": load_raccess_window_memory(),
    }

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # Panel A: runtime
    ax = axes[0]
    for method, (lengths, values) in runtime.items():
        if not lengths:
            continue
        lengths_arr = np.array(lengths)
        values_arr = np.array(values)
        ax.scatter(lengths_arr, values_arr, alpha=0.4, s=12, c=COLORS[method],
                   marker=MARKERS[method], edgecolors="none",
                   label=f"{method} (n={len(lengths_arr):,})")
        if method != "Raccess Exact" and len(lengths_arr) > 1:
            slope, intercept = linregress_slope_intercept(lengths_arr, values_arr)
            x_fit = np.linspace(lengths_arr.min(), lengths_arr.max(), 100)
            ax.plot(x_fit, slope * x_fit + intercept, color=COLORS[method], linewidth=1.8, alpha=0.9)
    ax.set_xlabel("Sequence length (nt)", fontsize=16)
    ax.set_ylabel("Runtime (seconds)", fontsize=16)
    ax.set_title("(A)", fontsize=18, fontweight="bold", loc="left")
    ax.legend(loc="upper left", fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis="both", labelsize=13)

    # Panel B: memory
    ax = axes[1]
    for method, (lengths, values) in memory.items():
        if not lengths:
            continue
        lengths_arr = np.array(lengths)
        values_arr = np.array(values)
        ax.scatter(lengths_arr, values_arr, alpha=0.4, s=12, c=COLORS[method],
                   marker=MARKERS[method], edgecolors="none",
                   label=f"{method} (n={len(lengths_arr):,})")
        if len(lengths_arr) > 1:
            if method == "Raccess Exact":
                coeffs = np.polyfit(lengths_arr, values_arr, 2)
                fit = np.poly1d(coeffs)
                x_fit = np.linspace(lengths_arr.min(), lengths_arr.max(), 100)
                ax.plot(x_fit, fit(x_fit), color=COLORS[method], linewidth=1.8, alpha=0.9)
            else:
                slope, intercept = linregress_slope_intercept(lengths_arr, values_arr)
                x_fit = np.linspace(lengths_arr.min(), lengths_arr.max(), 100)
                ax.plot(x_fit, slope * x_fit + intercept, color=COLORS[method], linewidth=1.8, alpha=0.9)
    ax.set_xlabel("Sequence length (nt)", fontsize=16)
    ax.set_ylabel("Peak RSS memory (MB)", fontsize=16)
    ax.set_title("(B)", fontsize=18, fontweight="bold", loc="left")
    ax.legend(loc="upper left", fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis="both", labelsize=13)

    fig.tight_layout()
    fig.savefig(OUTPUT, bbox_inches="tight")
    fig.savefig(OUTPUT.with_suffix(".png"), dpi=300, bbox_inches="tight")
    print(f"Figure saved to: {OUTPUT}")
    print(f"Figure saved to: {OUTPUT.with_suffix('.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
