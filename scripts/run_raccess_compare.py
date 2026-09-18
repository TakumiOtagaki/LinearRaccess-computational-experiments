#!/usr/bin/env python3
"""Run and compare Raccess and LinearRaccess exact profiles.

This script is for local/provenance-tracked comparisons. Raccess is an
externally installed reference and its source tree is intentionally not
committed to this repository.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from benchmark_common import (
    BenchmarkError,
    verify_binary_build,
    load_json_object,
    hardware_summary,
    run_timed,
    sha256_file,
    write_tsv,
)

RT_OUTPUT = 0.61633008


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--raccess-binary", required=True, type=Path)
    parser.add_argument("--linracc-binary", required=True, type=Path)
    parser.add_argument("--access-len", default="1,3,7,10")
    parser.add_argument("--max-span", type=int, default=0,
                        help="Use the maximum input length when omitted.")
    parser.add_argument("--energy-thr", type=float, default=1000.0)
    parser.add_argument("--energy", default="turner1999",
                        choices=["turner1999", "turner2004"])
    parser.add_argument("--linracc-engine", default="raccess",
                        choices=["raccess", "lincapr"])
    parser.add_argument("--beam", type=int, default=0)
    parser.add_argument("--source-revision",
                        default=None, help="Optional expected revision; must match the software lock and binary")
    parser.add_argument("--raccess-repo", required=True, type=Path)
    parser.add_argument("--raccess-commit", required=True)
    parser.add_argument("--raccess-patch-note", default="")
    return parser.parse_args()


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            name = line[1:].strip()
            chunks = []
        elif name is None:
            raise BenchmarkError(f"sequence occurs before header in {path}")
        else:
            chunks.append("".join(line.split()).upper().replace("T", "U"))
    if name is not None and chunks:
        records.append((name, "".join(chunks)))
    if not records:
        raise BenchmarkError(f"no FASTA records found in {path}")
    for name, sequence in records:
        if not sequence:
            raise BenchmarkError(f"empty sequence for {name!r}")
        invalid = next((base for base in sequence if base not in "ACGU"), None)
        if invalid is not None:
            raise BenchmarkError(f"unsupported base {invalid!r} in {name!r}")
    return records


def parse_raccess(path: Path) -> dict[str, dict[int, dict[int, float]]]:
    parsed: dict[str, dict[int, dict[int, float]]] = {}
    current: str | None = None
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            current = line[1:].strip()
            parsed.setdefault(current, {})
            continue
        if current is None or "\t" not in line:
            continue
        pos_text, payload = line.split("\t", 1)
        pos = int(pos_text)
        values: dict[int, float] = {}
        for chunk in payload.split(";"):
            chunk = chunk.strip()
            if not chunk:
                continue
            length_text, energy_text = chunk.split(",", 1)
            values[int(length_text)] = float(energy_text)
        parsed[current][pos] = values
    return parsed


def parse_linracc(path: Path) -> dict[str, dict[int, dict[int, float]]]:
    parsed: dict[str, dict[int, dict[int, float]]] = {}
    current: str | None = None
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.rstrip("\n")
        if not line:
            continue
        if line.startswith("# "):
            continue
        if line.startswith(">"):
            current = line[1:].strip()
            parsed.setdefault(current, {})
            continue
        if current is None or "\t" not in line:
            continue
        pos_text, payload = line.split("\t", 1)
        pos = int(pos_text)
        values: dict[int, float] = {}
        for chunk in payload.split(";"):
            chunk = chunk.strip()
            if not chunk:
                continue
            length_text, probability_text = chunk.split(",", 1)
            values[int(length_text)] = float(probability_text)
        parsed[current][pos] = values
    return parsed


def energy_from_probability(probability: float, rt: float = RT_OUTPUT) -> float:
    if probability <= 0.0:
        return math.inf
    return -rt * math.log(probability)


def probability_from_energy(energy: float, rt: float = RT_OUTPUT) -> float:
    if not math.isfinite(energy):
        return 0.0
    return math.exp(-energy / rt)


def percentile(values: list[float], percent: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    index = max(0, min(len(ordered) - 1, math.ceil(percent / 100.0 * len(ordered)) - 1))
    return ordered[index]


def summarize_values(values: list[float]) -> dict[str, float | None]:
    if not values:
        return {
            "count": 0,
            "min": None,
            "max": None,
            "mean": None,
            "median": None,
            "rms": None,
            "p95": None,
        }
    return {
        "count": len(values),
        "min": min(values),
        "max": max(values),
        "mean": sum(values) / len(values),
        "median": statistics.median(values),
        "rms": math.sqrt(sum(value * value for value in values) / len(values)),
        "p95": percentile(values, 95.0),
    }


def rank_values(values: list[float]) -> list[float]:
    indexed = sorted(range(len(values)), key=lambda index: values[index])
    ranks = [0.0] * len(values)
    index = 0
    while index < len(indexed):
        end = index
        while (
            end + 1 < len(indexed)
            and values[indexed[end + 1]] == values[indexed[index]]
        ):
            end += 1
        average_rank = 0.5 * (index + end) + 1.0
        for offset in range(index, end + 1):
            ranks[indexed[offset]] = average_rank
        index = end + 1
    return ranks


def spearman(x: list[float], y: list[float]) -> float | None:
    if len(x) != len(y) or len(x) < 2:
        return None
    x_rank = rank_values(x)
    y_rank = rank_values(y)
    x_mean = sum(x_rank) / len(x_rank)
    y_mean = sum(y_rank) / len(y_rank)
    numerator = sum((a - x_mean) * (b - y_mean) for a, b in zip(x_rank, y_rank))
    denominator = math.sqrt(
        sum((a - x_mean) ** 2 for a in x_rank)
        * sum((b - y_mean) ** 2 for b in y_rank)
    )
    if denominator == 0.0:
        return None
    return numerator / denominator


def compare_record(
    name: str,
    length: int,
    access_lens: list[int],
    raccess: dict[int, dict[int, float]],
    linracc: dict[int, dict[int, float]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    energy_deltas: list[float] = []
    probability_deltas: list[float] = []
    raccess_energies: list[float] = []
    linracc_energies: list[float] = []
    missing_raccess = 0
    missing_linracc = 0
    nonfinite = 0

    for pos in range(length):
        raccess_values = raccess.get(pos, {})
        linracc_values = linracc.get(pos, {})
        for access_len in access_lens:
            if pos + access_len > length:
                continue
            raccess_energy = raccess_values.get(access_len)
            linracc_probability = linracc_values.get(access_len)
            if raccess_energy is None:
                missing_raccess += 1
                continue
            if linracc_probability is None:
                missing_linracc += 1
                continue
            raccess_probability = probability_from_energy(raccess_energy)
            linracc_energy = energy_from_probability(linracc_probability)
            if not (math.isfinite(raccess_energy) and math.isfinite(linracc_energy)):
                nonfinite += 1
                rows.append(
                    {
                        "sequence": name,
                        "position": pos,
                        "access_len": access_len,
                        "raccess_energy": raccess_energy,
                        "raccess_probability": raccess_probability,
                        "linracc_probability": linracc_probability,
                        "linracc_energy": linracc_energy,
                        "energy_diff": math.inf,
                        "probability_diff": math.inf,
                    }
                )
                continue
            energy_diff = abs(linracc_energy - raccess_energy)
            probability_diff = abs(linracc_probability - raccess_probability)
            energy_deltas.append(energy_diff)
            probability_deltas.append(probability_diff)
            raccess_energies.append(raccess_energy)
            linracc_energies.append(linracc_energy)
            rows.append(
                {
                    "sequence": name,
                    "position": pos,
                    "access_len": access_len,
                    "raccess_energy": raccess_energy,
                    "raccess_probability": raccess_probability,
                    "linracc_probability": linracc_probability,
                    "linracc_energy": linracc_energy,
                    "energy_diff": energy_diff,
                    "probability_diff": probability_diff,
                }
            )

    summary = {
        "sequence": name,
        "length": length,
        "expected_cells": sum(
            max(0, length - access_len + 1) for access_len in access_lens
        ),
        "compared_cells": len(energy_deltas),
        "missing_raccess_cells": missing_raccess,
        "missing_linracc_cells": missing_linracc,
        "nonfinite_cells": nonfinite,
        "energy_diff": summarize_values(energy_deltas),
        "probability_diff": summarize_values(probability_deltas),
        "spearman_energy": spearman(raccess_energies, linracc_energies),
    }
    return rows, summary


def git_sha(path: Path) -> str:
    import subprocess

    completed = subprocess.run(
        ["git", "-C", str(path), "rev-parse", "HEAD"],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        return "unknown"
    return completed.stdout.strip()


def main() -> int:
    args = parse_args()
    if args.bundle.exists():
        raise BenchmarkError(f"result bundle already exists: {args.bundle}")

    expected = load_json_object(SCRIPT_DIR.parents[1] / "software.lock.json")["software"]["commit"]
    if args.source_revision is not None and args.source_revision != expected:
        raise BenchmarkError("requested source revision differs from software lock")
    build_info = verify_binary_build(args.linracc_binary.resolve(), expected)
    args.source_revision = build_info["source_revision"]
    records = parse_fasta(args.input)
    max_length = max(len(sequence) for _, sequence in records)
    max_span = args.max_span or max_length
    if max_span < max_length:
        raise BenchmarkError(
            f"max_span {max_span} is smaller than the longest input {max_length}"
        )
    access_lens = [int(value) for value in args.access_len.split(",") if value]
    if not access_lens or any(value <= 0 for value in access_lens):
        raise BenchmarkError("access lengths must be positive integers")
    if len(set(access_lens)) != len(access_lens):
        raise BenchmarkError("access lengths must be unique")
    args.bundle.mkdir(parents=True, exist_ok=False)
    (args.bundle / "raccess").mkdir()
    (args.bundle / "linearraccess").mkdir()
    (args.bundle / "analysis").mkdir()

    raccess_out = args.bundle / "raccess" / "raccess_out.txt"
    linracc_out = args.bundle / "linearraccess" / "linracc_out.txt"

    raccess_cmd = [
        str(args.raccess_binary),
        f"-outfile={raccess_out}",
        f"-seqfile={args.input}",
        f"-access_len={args.access_len}",
        f"-max_span={max_span}",
        f"-energy_thr={args.energy_thr}",
    ]
    linracc_cmd = [
        str(args.linracc_binary),
        f"-seqfile={args.input}",
        f"-outfile={linracc_out}",
        f"-access_len={args.access_len}",
        f"-engine={args.linracc_engine}",
        f"-energy={args.energy}",
        f"-beam={args.beam}",
        f"-c_multi={max_span}",
        f"-c_hairpin={max_span}",
        "-no-fast-logsumexp",
        "-no-normalize",
        "-probabilities",
        "-metadata",
    ]

    raccess_elapsed, raccess_peak_rss_mb = run_timed(raccess_cmd)
    linracc_elapsed, linracc_peak_rss_mb = run_timed(linracc_cmd)

    raccess_parsed = parse_raccess(raccess_out)
    linracc_parsed = parse_linracc(linracc_out)

    all_rows: list[dict[str, Any]] = []
    per_sequence: list[dict[str, Any]] = []
    all_energy_deltas: list[float] = []
    all_probability_deltas: list[float] = []
    all_raccess_energies: list[float] = []
    all_linracc_energies: list[float] = []

    for name, sequence in records:
        rows, summary = compare_record(
            name,
            len(sequence),
            access_lens,
            raccess_parsed.get(name, {}),
            linracc_parsed.get(name, {}),
        )
        all_rows.extend(rows)
        per_sequence.append(summary)
        for row in rows:
            if math.isfinite(row["energy_diff"]):
                all_energy_deltas.append(row["energy_diff"])
                all_probability_deltas.append(row["probability_diff"])
                all_raccess_energies.append(row["raccess_energy"])
                all_linracc_energies.append(row["linracc_energy"])

    if all_rows:
        fieldnames = [
            "sequence",
            "position",
            "access_len",
            "raccess_energy",
            "raccess_probability",
            "linracc_probability",
            "linracc_energy",
            "energy_diff",
            "probability_diff",
        ]
        write_tsv(args.bundle / "analysis" / "comparison.tsv", all_rows)
    else:
        (args.bundle / "analysis" / "comparison.tsv").write_text(
            "sequence\tposition\taccess_len\traccess_energy\traccess_probability\tlinracc_probability\tlinracc_energy\tenergy_diff\tprobability_diff\n"
        )

    summary = {
        "access_lens": access_lens,
        "energy": args.energy,
        "beam": args.beam,
        "max_span": max_span,
        "energy_thr": args.energy_thr,
        "total_rows": len(all_rows),
        "compared_cells": len(all_energy_deltas),
        "energy_diff": summarize_values(all_energy_deltas),
        "probability_diff": summarize_values(all_probability_deltas),
        "spearman_energy": spearman(all_raccess_energies, all_linracc_energies),
        "per_sequence": per_sequence,
        "runtime": {
            "raccess_seconds": raccess_elapsed,
            "raccess_peak_rss_mb": raccess_peak_rss_mb,
            "linracc_seconds": linracc_elapsed,
            "linracc_peak_rss_mb": linracc_peak_rss_mb,
        },
    }
    summary["comparison_complete"] = all(
        item["expected_cells"] > 0
        and item["compared_cells"] == item["expected_cells"]
        for item in per_sequence
    )
    (args.bundle / "analysis" / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (args.bundle / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    provenance = {
        "schema_version": 1,
        "bundle": str(args.bundle),
        "input": {
            "path": str(args.input),
            "sha256": sha256_file(args.input),
            "records": [
                {"id": name, "length": len(sequence)} for name, sequence in records
            ],
        },
        "access_lens": access_lens,
        "settings": {
            "energy": args.energy,
            "beam": args.beam,
            "max_span": max_span,
            "energy_thr": args.energy_thr,
            "linracc_engine": args.linracc_engine,
            "linracc_normalize": False,
            "linracc_logsumexp": "exact",
            "raccess_definition": "max_span>=N, full-span exact Raccess reference",
        },
        "raccess": {
            "binary": str(args.raccess_binary),
            "binary_sha256": sha256_file(args.raccess_binary),
            "repo": str(args.raccess_repo),
            "commit": args.raccess_commit,
            "checked_out_commit": git_sha(args.raccess_repo),
            "patch_note": args.raccess_patch_note,
            "command": " ".join(raccess_cmd),
        },
        "linracc": {
            "binary": str(args.linracc_binary),
            "binary_sha256": sha256_file(args.linracc_binary),
            "source_revision": args.source_revision,
            "build_info": build_info,
            "command": " ".join(linracc_cmd),
        },
        "outputs": {
            "raccess": str(raccess_out),
            "raccess_sha256": sha256_file(raccess_out),
            "linracc": str(linracc_out),
            "linracc_sha256": sha256_file(linracc_out),
        },
        "hardware": hardware_summary(),
    }
    (args.bundle / "provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(json.dumps(summary, indent=2, sort_keys=True))
    if not summary["comparison_complete"]:
        print("error: incomplete or nonfinite comparison; inspect per_sequence diagnostics", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
