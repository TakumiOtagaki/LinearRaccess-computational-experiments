#!/usr/bin/env python3
"""Compute the aggregate summary.json for benchmark/results/bprna_accuracy_comparison,
matching Table 1's methodology exactly (12,872 bpRNA sequences, 100-1,000 nt,
window lengths w in {1,3,7,10}, median absolute error per length bin).

This bundle previously tracked only per-sequence run metadata (*_summary.json:
status/memory/elapsed, no accessibility values) with the raw per-cell outputs
(exact/, linear_b100/, window_w100/) gitignored -- so a fresh clone had no
machine-readable aggregate supporting Table 1's numbers. This script produces
that aggregate, using the same binning and error definition as
benchmark/scripts/rnaplfold_compare.py (which already covers Table 1's
RNAplfold column).

Usage:
    python3 benchmark/scripts/compute_bprna_accuracy_summary.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "benchmark/scripts"))

from run_raccess_compare import parse_raccess, parse_linracc, probability_from_energy

RESULTS_DIR = REPO_ROOT / "benchmark/results/bprna_accuracy_comparison"
EXACT_DIR = RESULTS_DIR / "exact"
B100_DIR = RESULTS_DIR / "linear_b100"
WINDOW_DIR = RESULTS_DIR / "window_w100"
SEQ_LIST = REPO_ROOT / "benchmark/data/bprna_100_1000_list.txt"
ACCESS_LENS = (1, 3, 7, 10)

BINS = [(100, 300), (300, 500), (500, 700), (700, 1001)]
BIN_LABELS = ["100-300", "300-500", "500-700", "700-1000"]


def bin_for_length(length: int) -> int | None:
    for idx, (lo, hi) in enumerate(BINS):
        if lo <= length < hi:
            return idx
    return None


def all_errors(exact_map: dict, approx_map: dict, sequence_length: int) -> list[float]:
    errors = []
    for length in ACCESS_LENS:
        for pos in range(sequence_length - length + 1):
            try:
                exact_energy = exact_map[pos][length]
                p_approx = approx_map[pos][length]
            except KeyError as exc:
                raise ValueError(f"Missing cell at start={pos}, length={length}") from exc
            p_exact = probability_from_energy(exact_energy)
            if not math.isfinite(p_exact) or not math.isfinite(p_approx):
                raise ValueError(f"Non-finite probability at start={pos}, length={length}")
            errors.append(abs(p_exact - p_approx))
    return errors


def median(values: list[float]) -> float | None:
    if not values:
        return None
    s = sorted(values)
    n = len(s)
    mid = n // 2
    if n % 2 == 0:
        return (s[mid - 1] + s[mid]) / 2
    return s[mid]


def summarize(method_name: str, approx_dir: Path, is_linracc: bool, seqs: list[tuple[str, int]]) -> dict:
    errors_by_bin: list[list[float]] = [[] for _ in BINS]
    n_by_bin = [0] * len(BINS)
    n_missing = 0
    n_total_cells = 0

    for seq_id, length in seqs:
        b = bin_for_length(length)
        if b is None:
            continue
        exact_path = EXACT_DIR / f"{seq_id}.txt"
        approx_path = approx_dir / f"{seq_id}.txt"
        if not exact_path.exists() or not approx_path.exists():
            n_missing += 1
            continue
        n_by_bin[b] += 1
        exact_map = parse_raccess(exact_path).get(seq_id, {})
        if is_linracc:
            approx_map = parse_linracc(approx_path).get(seq_id, {})
        else:
            raw_map = parse_raccess(approx_path).get(seq_id, {})
            approx_map = {
                pos: {length: probability_from_energy(e) for length, e in lens.items()}
                for pos, lens in raw_map.items()
            }
        errs = all_errors(exact_map, approx_map, length)
        errors_by_bin[b].extend(errs)
        n_total_cells += len(errs)

    result = {
        "method": method_name,
        "n_missing_sequences": n_missing,
        "n_cells_compared": n_total_cells,
        "by_bin": {},
    }
    pooled: list[float] = []
    for idx, label in enumerate(BIN_LABELS):
        result["by_bin"][label] = {
            "n_sequences": n_by_bin[idx],
            "n_cells": len(errors_by_bin[idx]),
            "median_absolute_error": median(errors_by_bin[idx]),
        }
        pooled.extend(errors_by_bin[idx])
    result["pooled"] = {
        "n_sequences": sum(n_by_bin),
        "n_cells": len(pooled),
        "median_absolute_error": median(pooled),
    }
    return result


def main() -> int:
    seqs = []
    for line in SEQ_LIST.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        seq_id, length = line.split()
        seqs.append((seq_id, int(length)))

    linracc_summary = summarize("LinearRaccess (beam=100)", B100_DIR, True, seqs)
    window_summary = summarize("Raccess Window (W=100)", WINDOW_DIR, False, seqs)

    output = {
        "schema_version": 1,
        "description": (
            "Aggregate accuracy statistics for Table 1 / Figure 2: LinearRaccess "
            "(beam=100) and Raccess Window (W=100) median absolute error against "
            "Raccess Exact, by sequence-length bin, on 12,872 bpRNA-1m(90) "
            "sequences (100-1,000 nt) at window lengths w in {1,3,7,10}. "
            "Computed from the raw per-cell outputs in exact/, linear_b100/, "
            "and window_w100/ (gitignored; regenerable from software.lock.json's "
            "pinned LinearRaccess revision and the Raccess reference commit "
            "recorded below)."
        ),
        "software": {
            "linearraccess_source_revision": "892c93c80f4081b5a30f27901e7c16c25cac2088",
            "raccess_source_revision": "ad29d098b62606aa4051951469780b4e7be54536",
        },
        "bin_labels": BIN_LABELS,
        "access_lens": list(ACCESS_LENS),
        "sequence_list": "benchmark/data/bprna_100_1000_list.txt",
        "results": {
            "linearraccess_b100": linracc_summary,
            "raccess_window_w100": window_summary,
        },
    }

    out_path = RESULTS_DIR / "summary.json"
    out_path.write_text(json.dumps(output, indent=2, sort_keys=False) + "\n")
    print(f"Wrote {out_path}")
    print(f"LinearRaccess pooled median: {linracc_summary['pooled']['median_absolute_error']:.6g}")
    print(f"Raccess Window pooled median: {window_summary['pooled']['median_absolute_error']:.6g}")
    for label in BIN_LABELS:
        print(f"  [{label}] LinearRaccess={linracc_summary['by_bin'][label]['median_absolute_error']:.6g} "
              f"Window={window_summary['by_bin'][label]['median_absolute_error']:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
