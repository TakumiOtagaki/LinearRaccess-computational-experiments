#!/usr/bin/env python3
"""Compare RNAplfold's -u accessibility output against the Raccess Exact
reference, matching Table 1's methodology (12,872 bpRNA sequences, 100-1,000
nt, window lengths w in {1,3,7,10}).

RNAplfold computes the same conceptual quantity as Raccess (the probability
that a window is entirely unpaired) via a fundamentally different algorithm:
a sliding local-folding window of size -W, averaged over all windows
containing the region, rather than Raccess's single global ensemble with a
span cutoff. Settings here (-W 100 -L 100 -u 10) match Raccess Window's
W=100 span limit and the paper's tested window lengths.

Usage:
    python3 benchmark/scripts/rnaplfold_compare.py

Reads benchmark/results/rnaplfold_turner1999_w100_u10/lunp/*_lunp and
benchmark/results/bprna_accuracy_comparison/exact/*.txt, and writes
summary.json in the results directory. Use --results-dir for a different
bundle; its provenance.json must identify the settings used.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "benchmark/scripts"))

from run_raccess_compare import parse_raccess, probability_from_energy

RESULTS_DIR = REPO_ROOT / "benchmark/results/rnaplfold_turner1999_w100_u10"
EXACT_DIR = REPO_ROOT / "benchmark/results/bprna_accuracy_comparison/exact"
SEQ_LIST = REPO_ROOT / "benchmark/data/bprna_100_1000_list.txt"
ACCESS_LENS = (1, 3, 7, 10)

BINS = [(100, 300), (300, 500), (500, 700), (700, 1001)]
BIN_LABELS = ["100-300", "300-500", "500-700", "700-1000"]


def parse_lunp(path: Path) -> dict[int, dict[int, float]]:
    """Parse an RNAplfold _lunp file into {start_0idx: {length: prob}}.

    Row i (1-indexed), column l: probability that positions [i-l+1, i]
    (1-indexed, inclusive) are all unpaired. Converted to this repo's
    0-indexed start-anchored convention: start_0idx = i - l, matching
    parse_raccess/parse_linracc's {pos: {length: value}} shape.
    """
    result: dict[int, dict[int, float]] = {}
    lines = path.read_text().splitlines()
    for line in lines:
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        try:
            i = int(parts[0])
        except ValueError:
            continue
        for l_idx, cell in enumerate(parts[1:], start=1):
            cell = cell.strip()
            if cell == "NA" or cell == "":
                continue
            try:
                prob = float(cell)
            except ValueError:
                continue
            start_0idx = i - l_idx
            if start_0idx < 0:
                continue
            result.setdefault(start_0idx, {})[l_idx] = prob
    return result


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
            if not 0 <= p_approx <= 1:
                raise ValueError(f"RNAplfold probability outside [0,1]: {p_approx}")
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    args = parser.parse_args()
    provenance = json.loads((args.results_dir / "provenance.json").read_text())
    lunp_dir = args.results_dir / "lunp"
    seqs = []
    for line in SEQ_LIST.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        seq_id, length = line.split()
        seqs.append((seq_id, int(length)))

    errors_by_bin: list[list[float]] = [[] for _ in BINS]
    n_by_bin = [0] * len(BINS)
    n_missing = 0
    n_total_cells = 0

    for seq_id, length in seqs:
        b = bin_for_length(length)
        if b is None:
            continue
        lunp_path = lunp_dir / f"{seq_id}_lunp"
        exact_path = EXACT_DIR / f"{seq_id}.txt"
        if not lunp_path.exists() or not exact_path.exists():
            n_missing += 1
            continue
        n_by_bin[b] += 1
        rnaplfold_map = parse_lunp(lunp_path)
        exact_map = parse_raccess(exact_path).get(seq_id, {})
        try:
            errs = all_errors(exact_map, rnaplfold_map, length)
        except ValueError as exc:
            raise ValueError(f"{seq_id}: {exc}") from exc
        errors_by_bin[b].extend(errs)
        n_total_cells += len(errs)

    print(f"n_missing sequences: {n_missing}")
    print(f"n_cells_compared: {n_total_cells}")
    if n_missing:
        raise ValueError(f"Incomplete comparison: {n_missing} sequences missing")

    result = {
        "evidence_scope": (
            "RNAplfold (-W 100 -L 100 -u 10, ViennaRNA 2.7.2) accuracy "
            "comparison against Raccess Exact, matching Table 1's 12,872 "
            "bpRNA sequences (100-1,000 nt) and window lengths w in "
            "{1,3,7,10}. RNAplfold computes the same accessibility "
            "quantity via a different algorithm: local sliding-window "
            "folding averaged over all windows containing the region, "
            "rather than Raccess's single global ensemble with a span "
            "cutoff."
        ),
        "settings": provenance["settings"],
        "tool": provenance["tool"],
        "n_missing_sequences": n_missing,
        "n_cells_compared": n_total_cells,
        "by_bin": {},
    }

    all_errors_pooled: list[float] = []
    for idx, label in enumerate(BIN_LABELS):
        errs = errors_by_bin[idx]
        all_errors_pooled.extend(errs)
        result["by_bin"][label] = {
            "n_sequences": n_by_bin[idx],
            "n_cells": len(errs),
            "median_absolute_error": median(errs),
        }
        print(f"{label}: n={n_by_bin[idx]}, median_error={median(errs):.3e}" if errs else f"{label}: n=0")

    result["overall"] = {
        "n_sequences": sum(n_by_bin),
        "n_cells": len(all_errors_pooled),
        "median_absolute_error": median(all_errors_pooled),
    }
    print(f"Overall: n={sum(n_by_bin)}, median_error={median(all_errors_pooled):.3e}")

    output = args.results_dir / "summary.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
