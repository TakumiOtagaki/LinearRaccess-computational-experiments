#!/usr/bin/env python3
"""Extend benchmark/results/bprna_w_cap_sweep_300_1000/summary.json with the
p95/max relative-error statistics (Table S4) and the absolute-error
statistics (Supplementary S4 prose) that were previously computed only
ad hoc and not tracked, per an independent review's finding that the
tracked summary supported only Table S4's 21 medians, not its other 42
p95/max entries nor the newly-added absolute-error counts and maxima.

Usage:
    python3 benchmark/scripts/compute_w_cap_sweep_summary.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "benchmark/scripts"))

from run_raccess_compare import parse_raccess, parse_linracc, probability_from_energy

ROOT = REPO_ROOT / "benchmark/results/bprna_w_cap_sweep_300_1000"
SEQ_LIST = REPO_ROOT / "benchmark/data/bprna_300_1000_list.txt"
ACCESS_LENS = [1, 3, 7, 10, 20, 30, 50]


def load_seqs() -> list[tuple[str, int]]:
    seqs = []
    for line in SEQ_LIST.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        sid, length = line.split()
        seqs.append((sid, int(length)))
    return seqs


def load_approx(method_dir: str, is_linracc: bool, sid: str) -> dict:
    path = ROOT / method_dir / f"{sid}.txt"
    if not path.exists():
        return {}
    if is_linracc:
        return parse_linracc(path).get(sid, {})
    raw = parse_raccess(path).get(sid, {})
    return {pos: {l: probability_from_energy(e) for l, e in lens.items()} for pos, lens in raw.items()}


def percentile(values: list[float], p: float) -> float | None:
    if not values:
        return None
    s = sorted(values)
    n = len(s)
    k = (n - 1) * p / 100
    f, c = math.floor(k), math.ceil(k)
    if f == c:
        return s[int(k)]
    return s[f] + (s[c] - s[f]) * (k - f)


def median(values: list[float]) -> float | None:
    return percentile(values, 50)


def stats_for(w: int, method_dir: str, is_linracc: bool, seqs: list[tuple[str, int]]) -> dict:
    rel_errs: list[float] = []
    abs_errs: list[float] = []
    n_abs_gt_05 = 0
    for sid, length in seqs:
        exact_path = ROOT / "exact" / f"{sid}.txt"
        if not exact_path.exists():
            continue
        exact_map = parse_raccess(exact_path).get(sid, {})
        approx_map = load_approx(method_dir, is_linracc, sid)
        if not approx_map:
            continue
        for pos in range(length - w + 1):
            try:
                pe = probability_from_energy(exact_map[pos][w])
                pa = approx_map[pos][w]
            except KeyError:
                continue
            if not math.isfinite(pe) or not math.isfinite(pa):
                continue
            ae = abs(pe - pa)
            abs_errs.append(ae)
            if ae > 0.05:
                n_abs_gt_05 += 1
            if pe > 0:
                rel_errs.append(ae / pe)
    return dict(
        n_cells=len(rel_errs),
        median_relative_error=median(rel_errs),
        p95_relative_error=percentile(rel_errs, 95),
        max_relative_error=max(rel_errs) if rel_errs else None,
        max_absolute_error=max(abs_errs) if abs_errs else None,
        n_absolute_error_gt_0_05=n_abs_gt_05,
    )


def main() -> int:
    seqs = load_seqs()
    summary_path = ROOT / "summary.json"
    summary = json.loads(summary_path.read_text())

    detailed: dict[str, dict] = {}
    for w in ACCESS_LENS:
        detailed[str(w)] = dict(
            beam100_cap60=stats_for(w, "linear_b100_c60", True, seqs),
            beam200_cap60=stats_for(w, "linear_b200_c60", True, seqs),
            raccess_window_w100=stats_for(w, "window_w100", False, seqs),
        )
        print(f"w={w}: beam100 p95={detailed[str(w)]['beam100_cap60']['p95_relative_error']:.6g} "
              f"max_abs={detailed[str(w)]['beam100_cap60']['max_absolute_error']:.6g} "
              f"n_abs>0.05={detailed[str(w)]['beam100_cap60']['n_absolute_error_gt_0_05']}")

    summary["table_s4_full_statistics"] = {
        "description": (
            "Full per-w, per-method statistics backing Table S4 (median, p95, "
            "max relative error) and the S4 prose's absolute-error claims "
            "(max absolute error, count of cells with absolute error > 0.05), "
            "computed directly from the raw per-cell outputs in exact/, "
            "linear_b100_c60/, linear_b200_c60/, and window_w100/."
        ),
        "by_w": detailed,
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Updated {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
