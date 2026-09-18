#!/usr/bin/env python3
"""Characterize LinearRaccess's *relative* error at low true accessibility
(rare accessible windows), the regime the Discussion's determinism argument
is about. Table 1 only reports *absolute* error, which is dominated by
high-probability cells and says nothing about relative error where the true
value is small -- this was flagged as a gap during review (Discussion
claims determinism matters most for rare events, without reporting whether
LinearRaccess is itself accurate there in relative terms).

Usage:
    python3 benchmark/scripts/rare_event_relative_error.py

Reads the same frozen benchmark/results/bprna_accuracy_comparison/{exact,
linear_b100} evidence used for Table 1 (all 12,872 sequences, no
resampling), and writes benchmark/results/rare_event_relative_error/
{provenance,summary}.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "benchmark/scripts"))

from run_raccess_compare import parse_raccess, parse_linracc, probability_from_energy

EXACT_DIR = REPO_ROOT / "benchmark/results/bprna_accuracy_comparison/exact"
B100_DIR = REPO_ROOT / "benchmark/results/bprna_accuracy_comparison/linear_b100"
SEQ_LIST = REPO_ROOT / "benchmark/data/bprna_100_1000_list.txt"
RESULTS_DIR = REPO_ROOT / "benchmark/results/rare_event_relative_error"
OUTPUT = RESULTS_DIR / "summary.json"
PROVENANCE = RESULTS_DIR / "provenance.json"

LOW, HIGH = 1e-8, 1e-5


def median(values: list[float]) -> float:
    s = sorted(values)
    n = len(s)
    mid = n // 2
    return (s[mid - 1] + s[mid]) / 2 if n % 2 == 0 else s[mid]


def percentile(values: list[float], p: float) -> float:
    s = sorted(values)
    idx = min(len(s) - 1, max(0, int(p / 100.0 * len(s))))
    return s[idx]


def main() -> int:
    seqs = [line.split()[0] for line in SEQ_LIST.read_text().splitlines() if line.strip()]

    rel_errors: list[float] = []
    n_total_cells = 0

    for seq_id in seqs:
        fname = seq_id + ".txt"
        ep = EXACT_DIR / fname
        bp = B100_DIR / fname
        if not (ep.exists() and bp.exists()):
            continue
        exact_map = parse_raccess(ep).get(seq_id, {})
        approx_map = parse_linracc(bp).get(seq_id, {})
        for pos, lens in exact_map.items():
            arow = approx_map.get(pos)
            if not arow:
                continue
            for length, exact_energy in lens.items():
                if length not in arow:
                    continue
                p_exact = probability_from_energy(exact_energy)
                n_total_cells += 1
                if LOW <= p_exact <= HIGH:
                    p_approx = arow[length]
                    rel_errors.append(abs(p_approx - p_exact) / p_exact)

    n_rare = len(rel_errors)
    frac_over_1 = sum(1 for x in rel_errors if x > 1.0) / n_rare

    provenance = {
        "schema_version": 1,
        "description": (
            "Relative error |P_approx - P_exact| / P_exact for LinearRaccess "
            "(beam=100) restricted to cells where the true (Raccess Exact) "
            "accessibility falls in a rare-event range -- Table 1/Fig. 2 only "
            "report absolute error, which is dominated by high-probability "
            "cells and says nothing about relative error where the true "
            "value is small. Supports Results \\S3.2."
        ),
        "reused_data": {
            "exact_reference": "benchmark/results/bprna_accuracy_comparison/exact",
            "linear_b100": "benchmark/results/bprna_accuracy_comparison/linear_b100",
            "sequence_list": "benchmark/data/bprna_100_1000_list.txt",
        },
        "true_probability_range": [LOW, HIGH],
        "n_sequences": len(seqs),
    }
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PROVENANCE.write_text(json.dumps(provenance, indent=2))

    summary = {
        "n_total_cells": n_total_cells,
        "n_cells_in_range": n_rare,
        "fraction_of_all_cells_in_range": n_rare / n_total_cells,
        "median_relative_error": median(rel_errors),
        "p90_relative_error": percentile(rel_errors, 90),
        "p99_relative_error": percentile(rel_errors, 99),
        "max_relative_error": max(rel_errors),
        "fraction_relative_error_over_1": frac_over_1,
    }
    OUTPUT.write_text(json.dumps(summary, indent=2))

    print(json.dumps({**provenance, **summary}, indent=2))
    print(f"\nWrote {PROVENANCE}")
    print(f"Wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
