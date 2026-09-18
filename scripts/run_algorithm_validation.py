#!/usr/bin/env python3
"""Measure the LinearRaccess approximation envelope against an internal exact mode."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any

from benchmark_common import (
    BenchmarkError,
    run_timed,
    git_snapshot,
    hardware_summary,
    load_json_object,
    parse_output_metadata,
    sha256_file,
    sha256_text,
    verify_software_checkout,
    verify_binary_build,
    write_tsv,
)


PAPER_ROOT = Path(__file__).resolve().parents[2]
LOCK_PATH = PAPER_ROOT / "software.lock.json"


def read_multi_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise BenchmarkError(f"cannot read FASTA {path}: {error}") from error
    for raw_line in lines:
        line = raw_line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            name = line[1:].strip()
            chunks = []
            if not name:
                raise BenchmarkError(f"empty FASTA identifier in {path}")
        elif name is None:
            raise BenchmarkError(f"sequence occurs before FASTA identifier in {path}")
        else:
            chunks.append("".join(line.split()).upper().replace("T", "U"))
    if name is not None:
        records.append((name, "".join(chunks)))
    if not records:
        raise BenchmarkError(f"FASTA contains no records: {path}")
    if len({name for name, _ in records}) != len(records):
        raise BenchmarkError(f"duplicate FASTA identifier in {path}")
    for name, sequence in records:
        if not sequence or any(base not in "ACGU" for base in sequence):
            raise BenchmarkError(f"invalid RNA sequence for {name!r} in {path}")
    return records


def positive_int(value: Any, label: str, *, allow_zero: bool = False) -> int:
    lower = 0 if allow_zero else 1
    if isinstance(value, bool) or not isinstance(value, int) or value < lower:
        qualifier = "non-negative" if allow_zero else "positive"
        raise BenchmarkError(f"{label} must be a {qualifier} integer")
    return value


def validate_config(raw: Any, label: str) -> dict[str, Any]:
    if not isinstance(raw, dict):
        raise BenchmarkError(f"{label} must be an object")
    config = dict(raw)
    for field in ("id", "family"):
        if not isinstance(config.get(field), str) or not config[field]:
            raise BenchmarkError(f"{label}.{field} must be a non-empty string")
    if config.get("engine") != "lincapr":
        raise BenchmarkError(f"{label}.engine must be lincapr")
    if config.get("energy") not in {"turner1999", "turner2004"}:
        raise BenchmarkError(f"{label}.energy must name a bundled Turner model")
    positive_int(config.get("beam"), f"{label}.beam", allow_zero=True)
    for field in ("c_multi", "c_hairpin"):
        value = config.get(field)
        if value != "sequence_length":
            positive_int(value, f"{label}.{field}")
    for field in ("fast_logsumexp", "normalize_profiles"):
        if not isinstance(config.get(field), bool):
            raise BenchmarkError(f"{label}.{field} must be boolean")
    return config


def validate_manifest(path: Path) -> dict[str, Any]:
    manifest = load_json_object(path)
    if manifest.get("schema_version") != 1:
        raise BenchmarkError("algorithm-validation manifest must use schema_version 1")
    if not isinstance(manifest.get("experiment_id"), str) or not manifest["experiment_id"]:
        raise BenchmarkError("manifest experiment_id must be a non-empty string")
    source = manifest.get("source")
    if not isinstance(source, dict):
        raise BenchmarkError("manifest source must be an object")
    relative = Path(str(source.get("path", "")))
    if not str(relative) or relative.is_absolute():
        raise BenchmarkError("source.path must be relative to the manifest")
    source_path = (path.parent / relative).resolve()
    if sha256_file(source_path) != source.get("file_sha256"):
        raise BenchmarkError("source FASTA SHA-256 mismatch")
    records = read_multi_fasta(source_path)
    expected = source.get("records")
    if not isinstance(expected, list) or len(expected) != len(records):
        raise BenchmarkError("source.records does not match FASTA record count")
    for index, ((name, sequence), item) in enumerate(zip(records, expected)):
        if not isinstance(item, dict):
            raise BenchmarkError(f"source.records[{index}] must be an object")
        observed = {
            "id": name,
            "length": len(sequence),
            "sequence_sha256": sha256_text(sequence),
        }
        if item != observed:
            raise BenchmarkError(
                f"source.records[{index}] mismatch: expected {item}, observed {observed}"
            )
    access_lens = manifest.get("access_lens")
    if not isinstance(access_lens, list) or not access_lens:
        raise BenchmarkError("access_lens must be a non-empty array")
    access_lens = [positive_int(value, "access_lens item") for value in access_lens]
    if access_lens != sorted(set(access_lens)):
        raise BenchmarkError("access_lens must be sorted and unique")
    if access_lens[-1] > min(len(sequence) for _, sequence in records):
        raise BenchmarkError("access length exceeds the shortest sequence")
    reference = validate_config(manifest.get("reference"), "reference")
    configurations_raw = manifest.get("configurations")
    if not isinstance(configurations_raw, list) or not configurations_raw:
        raise BenchmarkError("configurations must be a non-empty array")
    configurations = [
        validate_config(value, f"configurations[{index}]")
        for index, value in enumerate(configurations_raw)
    ]
    ids = [reference["id"], *(config["id"] for config in configurations)]
    if len(ids) != len(set(ids)):
        raise BenchmarkError("configuration identifiers must be unique")
    relative_floor = manifest.get("relative_error_floor")
    tolerance = manifest.get("invariant_tolerance")
    for value, label in ((relative_floor, "relative_error_floor"), (tolerance, "invariant_tolerance")):
        if not isinstance(value, (int, float)) or isinstance(value, bool) or value <= 0:
            raise BenchmarkError(f"{label} must be positive")
    return {
        "manifest": manifest,
        "source_path": source_path,
        "records": records,
        "access_lens": access_lens,
        "reference": reference,
        "configurations": configurations,
        "relative_error_floor": float(relative_floor),
        "invariant_tolerance": float(tolerance),
    }


def resolve_config(config: dict[str, Any], length: int) -> dict[str, Any]:
    resolved = dict(config)
    for field in ("c_multi", "c_hairpin"):
        if resolved[field] == "sequence_length":
            resolved[field] = length
    return resolved


def command_for(
    executable: Path,
    fasta_path: Path,
    output_path: Path,
    access_lens: list[int],
    config: dict[str, Any],
    software_commit: str,
) -> list[str]:
    command = [
        str(executable),
        f"-seqfile={fasta_path}",
        f"-outfile={output_path}",
        f"-access_len={','.join(map(str, access_lens))}",
        "-engine=lincapr",
        f"-energy={config['energy']}",
        f"-beam={config['beam']}",
        f"-c_multi={config['c_multi']}",
        f"-c_hairpin={config['c_hairpin']}",
        "-probabilities",
        "-metadata",
    ]
    if not config["normalize_profiles"]:
        command.append("-no-normalize")
    if not config["fast_logsumexp"]:
        command.append("-no-fast-logsumexp")
    return command


def parse_probabilities(
    path: Path, length: int, access_lens: list[int]
) -> dict[tuple[int, int], float]:
    values: dict[tuple[int, int], float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith(("#", ">")):
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            raise BenchmarkError(f"malformed probability row in {path}: {line!r}")
        start = int(fields[0])
        for item in fields[1].split(";"):
            if not item:
                continue
            window_raw, separator, probability_raw = item.partition(",")
            if not separator:
                raise BenchmarkError(f"malformed probability item in {path}: {item!r}")
            key = (start, int(window_raw))
            if key in values:
                raise BenchmarkError(f"duplicate probability {key} in {path}")
            probability = float(probability_raw)
            if not math.isfinite(probability):
                raise BenchmarkError(f"non-finite probability {key} in {path}")
            values[key] = probability
    expected = {
        (start, window)
        for window in access_lens
        for start in range(length - window + 1)
    }
    if set(values) != expected:
        raise BenchmarkError(
            f"probability grid mismatch in {path}: expected {len(expected)}, got {len(values)}"
        )
    return values


def check_metadata(
    metadata: dict[str, str],
    access_lens: list[int],
    config: dict[str, Any],
    software_commit: str,
) -> None:
    expected = {
        "linearraccess_metadata": "1",
        "source_revision": software_commit,
        "engine": "lincapr",
        "effective_energy_model": str(config["energy"]),
        "beam": str(config["beam"]),
        "c_multi": str(config["c_multi"]),
        "c_hairpin": str(config["c_hairpin"]),
        "logsumexp": "fast" if config["fast_logsumexp"] else "exact",
        "normalize_profiles": "true" if config["normalize_profiles"] else "false",
        "output_mode": "probability",
        "access_lens": ",".join(map(str, access_lens)),
    }
    mismatch = {
        key: {"expected": value, "observed": metadata.get(key)}
        for key, value in expected.items()
        if metadata.get(key) != value
    }
    if mismatch:
        raise BenchmarkError(f"LinearRaccess metadata mismatch: {mismatch}")


def percentile(values: list[float], percentage: float) -> float:
    if not values:
        raise BenchmarkError("cannot calculate a percentile of an empty sample")
    ordered = sorted(values)
    position = (len(ordered) - 1) * percentage / 100.0
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    cursor = 0
    while cursor < len(order):
        end = cursor + 1
        while end < len(order) and values[order[end]] == values[order[cursor]]:
            end += 1
        rank = (cursor + 1 + end) / 2.0
        for index in order[cursor:end]:
            ranks[index] = rank
        cursor = end
    return ranks


def pearson(left: list[float], right: list[float]) -> float | None:
    if len(left) != len(right) or len(left) < 2:
        return None
    mean_left = statistics.fmean(left)
    mean_right = statistics.fmean(right)
    numerator = sum(
        (x - mean_left) * (y - mean_right) for x, y in zip(left, right)
    )
    denominator = math.sqrt(
        sum((x - mean_left) ** 2 for x in left)
        * sum((y - mean_right) ** 2 for y in right)
    )
    return numerator / denominator if denominator else None


def summarize_errors(
    observed: dict[tuple[int, int], float],
    reference: dict[tuple[int, int], float],
    relative_floor: float,
) -> dict[str, Any]:
    keys = sorted(reference)
    absolute = [abs(observed[key] - reference[key]) for key in keys]
    relative = [
        abs(observed[key] - reference[key]) / reference[key]
        for key in keys
        if reference[key] >= relative_floor
    ]
    return {
        "value_count": len(keys),
        "median_absolute_error": statistics.median(absolute),
        "p95_absolute_error": percentile(absolute, 95.0),
        "maximum_absolute_error": max(absolute),
        "rmse": math.sqrt(statistics.fmean(value * value for value in absolute)),
        "relative_error_floor": relative_floor,
        "relative_error_count": len(relative),
        "median_relative_error": statistics.median(relative) if relative else None,
        "p95_relative_error": percentile(relative, 95.0) if relative else None,
        "maximum_relative_error": max(relative) if relative else None,
        "spearman_rho": pearson(
            average_ranks([observed[key] for key in keys]),
            average_ranks([reference[key] for key in keys]),
        ),
    }


def invariant_counts(
    values: dict[tuple[int, int], float],
    length: int,
    access_lens: list[int],
    tolerance: float,
) -> dict[str, int]:
    bound_violations = sum(
        value < -tolerance or value > 1.0 + tolerance for value in values.values()
    )
    monotonicity_violations = 0
    for shorter, longer in zip(access_lens, access_lens[1:]):
        offset = longer - shorter
        for start in range(length - longer + 1):
            long_value = values[(start, longer)]
            if long_value > values[(start, shorter)] + tolerance:
                monotonicity_violations += 1
            if long_value > values[(start + offset, shorter)] + tolerance:
                monotonicity_violations += 1
    return {
        "probability_bound_violations": bound_violations,
        "nested_window_monotonicity_violations": monotonicity_violations,
    }


def run_configuration(
    executable: Path,
    sequence_path: Path,
    output_path: Path,
    access_lens: list[int],
    config: dict[str, Any],
    software_commit: str,
) -> tuple[dict[tuple[int, int], float], float, str, float | None]:
    command = command_for(
        executable, sequence_path, output_path, access_lens, config, software_commit
    )
    seconds, peak_rss_mb = run_timed(command)
    check_metadata(parse_output_metadata(output_path), access_lens, config, software_commit)
    length = sum(
        len(line.strip())
        for line in sequence_path.read_text(encoding="ascii").splitlines()
        if not line.startswith(">")
    )
    return (
        parse_probabilities(output_path, length, access_lens),
        seconds,
        sha256_file(output_path),
        peak_rss_mb,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--software-dir", type=Path, required=True)
    parser.add_argument("--linracc", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest_path = args.manifest.resolve()
        validated = validate_manifest(manifest_path)
        paper_snapshot = git_snapshot(PAPER_ROOT)
        if paper_snapshot["status"]:
            raise BenchmarkError("paper repository must be clean before a benchmark run")
        lock = load_json_object(LOCK_PATH)
        software_snapshot = verify_software_checkout(args.software_dir.resolve(), lock)
        executable = args.linracc.resolve()
        if not executable.is_file() or not executable.stat().st_mode & 0o111:
            raise BenchmarkError(f"LinearRaccess executable is unavailable: {executable}")
        build_info = verify_binary_build(executable, software_snapshot["commit"])
        output_dir = args.output_dir.resolve()
        if output_dir.exists():
            raise BenchmarkError(f"result bundle already exists: {output_dir}")
        output_dir.parent.mkdir(parents=True, exist_ok=True)

        manifest = validated["manifest"]
        access_lens = validated["access_lens"]
        software_commit = software_snapshot["commit"]
        raw_rows: list[dict[str, Any]] = []
        comparison_rows: list[dict[str, Any]] = []
        all_values: dict[str, list[float]] = {
            config["id"]: [] for config in validated["configurations"]
        }
        all_references: dict[str, list[float]] = {
            config["id"]: [] for config in validated["configurations"]
        }
        staging = Path(
            tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent)
        )
        try:
            with tempfile.TemporaryDirectory(prefix="linearraccess_algorithm_") as raw:
                temporary = Path(raw)
                for sequence_index, (name, sequence) in enumerate(validated["records"]):
                    fasta_path = temporary / f"sequence_{sequence_index}.fa"
                    fasta_path.write_text(f">{name}\n{sequence}\n", encoding="ascii")
                    reference_config = resolve_config(validated["reference"], len(sequence))
                    reference_path = temporary / f"sequence_{sequence_index}_reference.txt"
                    reference, reference_seconds, reference_hash, reference_rss = run_configuration(
                        executable,
                        fasta_path,
                        reference_path,
                        access_lens,
                        reference_config,
                        software_commit,
                    )
                    reference_invariants = invariant_counts(
                        reference,
                        len(sequence),
                        access_lens,
                        validated["invariant_tolerance"],
                    )
                    raw_rows.append(
                        {
                            "sequence_id": name,
                            "sequence_length": len(sequence),
                            "sequence_sha256": sha256_text(sequence),
                            "configuration": reference_config["id"],
                            "family": reference_config["family"],
                            "beam": reference_config["beam"],
                            "c_multi": reference_config["c_multi"],
                            "c_hairpin": reference_config["c_hairpin"],
                            "fast_logsumexp": int(reference_config["fast_logsumexp"]),
                            "seconds": reference_seconds,
                            "peak_rss_mb": reference_rss,
                            "output_sha256": reference_hash,
                            **reference_invariants,
                            "median_absolute_error": 0.0,
                            "p95_absolute_error": 0.0,
                            "maximum_absolute_error": 0.0,
                            "rmse": 0.0,
                            "relative_error_count": len(reference),
                            "median_relative_error": 0.0,
                            "p95_relative_error": 0.0,
                            "maximum_relative_error": 0.0,
                            "spearman_rho": 1.0,
                        }
                    )
                    for raw_config in validated["configurations"]:
                        config = resolve_config(raw_config, len(sequence))
                        result_path = temporary / (
                            f"sequence_{sequence_index}_{config['id']}.txt"
                        )
                        observed, seconds, output_hash, peak_rss_mb = run_configuration(
                            executable,
                            fasta_path,
                            result_path,
                            access_lens,
                            config,
                            software_commit,
                        )
                        metrics = summarize_errors(
                            observed, reference, validated["relative_error_floor"]
                        )
                        invariants = invariant_counts(
                            observed,
                            len(sequence),
                            access_lens,
                            validated["invariant_tolerance"],
                        )
                        keys = sorted(reference)
                        all_values[config["id"]].extend(observed[key] for key in keys)
                        all_references[config["id"]].extend(reference[key] for key in keys)
                        for start, window in keys:
                            reference_value = reference[(start, window)]
                            observed_value = observed[(start, window)]
                            absolute_error = abs(observed_value - reference_value)
                            comparison_rows.append(
                                {
                                    "sequence_id": name,
                                    "sequence_length": len(sequence),
                                    "configuration": config["id"],
                                    "family": config["family"],
                                    "start_0based": start,
                                    "window_length": window,
                                    "reference_probability": reference_value,
                                    "observed_probability": observed_value,
                                    "absolute_error": absolute_error,
                                    "relative_error": (
                                        absolute_error / reference_value
                                        if reference_value
                                        >= validated["relative_error_floor"]
                                        else None
                                    ),
                                }
                            )
                        raw_rows.append(
                            {
                                "sequence_id": name,
                                "sequence_length": len(sequence),
                                "sequence_sha256": sha256_text(sequence),
                                "configuration": config["id"],
                                "family": config["family"],
                                "beam": config["beam"],
                                "c_multi": config["c_multi"],
                                "c_hairpin": config["c_hairpin"],
                                "fast_logsumexp": int(config["fast_logsumexp"]),
                                "seconds": seconds,
                                "peak_rss_mb": peak_rss_mb,
                                "output_sha256": output_hash,
                                **invariants,
                                "median_absolute_error": metrics["median_absolute_error"],
                                "p95_absolute_error": metrics["p95_absolute_error"],
                                "maximum_absolute_error": metrics["maximum_absolute_error"],
                                "rmse": metrics["rmse"],
                                "relative_error_count": metrics["relative_error_count"],
                                "median_relative_error": metrics["median_relative_error"],
                                "p95_relative_error": metrics["p95_relative_error"],
                                "maximum_relative_error": metrics["maximum_relative_error"],
                                "spearman_rho": metrics["spearman_rho"],
                            }
                        )
                        print(
                            f"{name}\t{config['id']}\t"
                            f"max_abs={metrics['maximum_absolute_error']:.6g}\t"
                            f"{seconds:.3f}s",
                            flush=True,
                        )

            summary_rows: list[dict[str, Any]] = []
            for config in validated["configurations"]:
                config_id = config["id"]
                observed_values = all_values[config_id]
                reference_values = all_references[config_id]
                observed = dict(enumerate(observed_values))
                reference = dict(enumerate(reference_values))
                metrics = summarize_errors(
                    observed, reference, validated["relative_error_floor"]
                )
                matching = [row for row in raw_rows if row["configuration"] == config_id]
                summary_rows.append(
                    {
                        "configuration": config_id,
                        "family": config["family"],
                        "sequence_count": len(validated["records"]),
                        **metrics,
                        "probability_bound_violations": sum(
                            row["probability_bound_violations"] for row in matching
                        ),
                        "nested_window_monotonicity_violations": sum(
                            row["nested_window_monotonicity_violations"] for row in matching
                        ),
                        "total_seconds": sum(row["seconds"] for row in matching),
                        "max_peak_rss_mb": max((row["peak_rss_mb"] for row in matching if row["peak_rss_mb"] is not None), default=None),
                    }
                )

            write_tsv(staging / "raw.tsv", raw_rows)
            write_tsv(staging / "comparisons.tsv", comparison_rows)
            write_tsv(staging / "summary.tsv", summary_rows)
            (staging / "summary.json").write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "experiment_id": manifest["experiment_id"],
                        "reference": validated["reference"]["id"],
                        "evidence_scope": (
                            "internal same-model exact-mode comparison on manifest-declared RNA; "
                            "not an independent Raccess-oracle or biological validation"
                        ),
                        "measurements": summary_rows,
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n",
                encoding="utf-8",
            )
            shutil.copyfile(manifest_path, staging / "manifest.json")
            (staging / "provenance.json").write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                        "experiment_id": manifest["experiment_id"],
                        "paper_repository": paper_snapshot,
                        "software_repository": software_snapshot,
                        "software_lock_sha256": sha256_file(LOCK_PATH),
                        "manifest_path": str(manifest_path.relative_to(PAPER_ROOT)),
                        "manifest_sha256": sha256_file(manifest_path),
                        "source_path": str(validated["source_path"].relative_to(PAPER_ROOT)),
                        "source_file_sha256": sha256_file(validated["source_path"]),
                        "linracc_path": str(executable),
                        "linracc_sha256": sha256_file(executable),
                        "build_info": build_info,
                        "hardware": hardware_summary(),
                        "argv": sys.argv,
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n",
                encoding="utf-8",
            )
            staging.rename(output_dir)
        except Exception:
            shutil.rmtree(staging, ignore_errors=True)
            raise
    except (BenchmarkError, OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    print(output_dir / "summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
