#!/usr/bin/env python3
"""Run a pinned, fixed-parameter LinearRaccess scaling benchmark."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import shutil
import statistics
import sys
import tempfile
from pathlib import Path
from typing import Any

from benchmark_common import (
    BenchmarkError,
    fit_log_log,
    git_snapshot,
    hardware_summary,
    load_json_object,
    parse_output_metadata,
    read_single_fasta,
    run_timed,
    sha256_file,
    sha256_text,
    verify_software_checkout,
    verify_binary_build,
    write_tsv,
)


PAPER_ROOT = Path(__file__).resolve().parents[2]
LOCK_PATH = PAPER_ROOT / "software.lock.json"


def positive_int(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise BenchmarkError(f"{label} must be a positive integer")
    return value


def validate_manifest(path: Path, tier_name: str) -> dict[str, Any]:
    manifest = load_json_object(path)
    if manifest.get("schema_version") != 1:
        raise BenchmarkError("scaling manifest must use schema_version 1")
    experiment_id = manifest.get("experiment_id")
    if not isinstance(experiment_id, str) or not experiment_id:
        raise BenchmarkError("manifest experiment_id must be a non-empty string")

    source = manifest.get("source")
    if not isinstance(source, dict):
        raise BenchmarkError("manifest source must be an object")
    source_relative = Path(str(source.get("path", "")))
    if not str(source_relative) or source_relative.is_absolute():
        raise BenchmarkError("source.path must be a non-empty path relative to the manifest")
    source_path = (path.parent / source_relative).resolve()
    if not source_path.is_file():
        raise BenchmarkError(f"source FASTA does not exist: {source_path}")
    name, sequence = read_single_fasta(source_path)
    expected_length = positive_int(source.get("length"), "source.length")
    if len(sequence) != expected_length:
        raise BenchmarkError(
            f"source length mismatch: expected {expected_length}, got {len(sequence)}"
        )
    actual_file_sha = sha256_file(source_path)
    if source.get("file_sha256") != actual_file_sha:
        raise BenchmarkError(
            f"source file SHA-256 mismatch: expected {source.get('file_sha256')}, "
            f"got {actual_file_sha}"
        )
    actual_sequence_sha = sha256_text(sequence)
    if source.get("sequence_sha256") != actual_sequence_sha:
        raise BenchmarkError(
            "normalized source sequence SHA-256 mismatch: "
            f"expected {source.get('sequence_sha256')}, got {actual_sequence_sha}"
        )

    access_lens = manifest.get("access_lens")
    if not isinstance(access_lens, list) or not access_lens:
        raise BenchmarkError("manifest access_lens must be a non-empty array")
    normalized_lens = [positive_int(value, "access_lens item") for value in access_lens]
    if normalized_lens != sorted(set(normalized_lens)):
        raise BenchmarkError("access_lens must be unique and sorted")

    configuration = manifest.get("configuration")
    if not isinstance(configuration, dict):
        raise BenchmarkError("manifest configuration must be an object")
    if configuration.get("engine") not in {"lincapr", "raccess"}:
        raise BenchmarkError("configuration.engine must be lincapr or raccess")
    if configuration.get("engine") == "lincapr" and configuration.get("energy") not in {
        "turner1999",
        "turner2004",
    }:
        raise BenchmarkError("the lincapr engine requires a named Turner energy model")
    for field in ("beam", "c_multi", "c_hairpin"):
        positive_int(configuration.get(field), f"configuration.{field}")
    for field in ("fast_logsumexp", "normalize_profiles"):
        if not isinstance(configuration.get(field), bool):
            raise BenchmarkError(f"configuration.{field} must be boolean")

    tiers = manifest.get("tiers")
    if not isinstance(tiers, dict) or tier_name not in tiers:
        available = ", ".join(sorted(tiers)) if isinstance(tiers, dict) else "none"
        raise BenchmarkError(f"unknown tier {tier_name!r}; available: {available}")
    tier = tiers[tier_name]
    if not isinstance(tier, dict):
        raise BenchmarkError(f"tier {tier_name!r} must be an object")
    lengths = tier.get("lengths")
    if not isinstance(lengths, list) or len(lengths) < 2:
        raise BenchmarkError("a scaling tier requires at least two lengths")
    normalized_lengths = [positive_int(value, "tier length") for value in lengths]
    if normalized_lengths != sorted(set(normalized_lengths)):
        raise BenchmarkError("tier lengths must be unique and sorted")
    if normalized_lengths[-1] > len(sequence):
        raise BenchmarkError("tier length exceeds the normalized source length")
    if normalized_lengths[0] < normalized_lens[-1]:
        raise BenchmarkError("every prefix must be at least the maximum access length")
    warmups = tier.get("warmup_runs")
    if isinstance(warmups, bool) or not isinstance(warmups, int) or warmups < 0:
        raise BenchmarkError("tier warmup_runs must be a non-negative integer")
    positive_int(tier.get("repetitions"), "tier.repetitions")

    return {
        "manifest": manifest,
        "source_path": source_path,
        "source_name": name,
        "sequence": sequence,
        "lengths": normalized_lengths,
        "access_lens": normalized_lens,
        "configuration": configuration,
        "tier": tier,
    }


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
        f"-engine={config['engine']}",
        f"-beam={config['beam']}",
        f"-c_multi={config['c_multi']}",
        f"-c_hairpin={config['c_hairpin']}",
        "-probabilities",
        "-metadata",
    ]
    if config["engine"] == "lincapr":
        command.append(f"-energy={config['energy']}")
    if not config["normalize_profiles"]:
        command.append("-no-normalize")
    if not config["fast_logsumexp"]:
        command.append("-no-fast-logsumexp")
    return command


def check_metadata(
    metadata: dict[str, str],
    access_lens: list[int],
    config: dict[str, Any],
    software_commit: str,
) -> None:
    expected = {
        "linearraccess_metadata": "1",
        "source_revision": software_commit,
        "engine": str(config["engine"]),
        "beam": str(config["beam"]),
        "c_multi": str(config["c_multi"]),
        "c_hairpin": str(config["c_hairpin"]),
        "logsumexp": "fast" if config["fast_logsumexp"] else "exact",
        "normalize_profiles": "true" if config["normalize_profiles"] else "false",
        "output_mode": "probability",
        "access_lens": ",".join(map(str, access_lens)),
    }
    if config["engine"] == "lincapr":
        expected["effective_energy_model"] = str(config["energy"])
    mismatches = {
        key: {"expected": value, "observed": metadata.get(key)}
        for key, value in expected.items()
        if metadata.get(key) != value
    }
    if mismatches:
        raise BenchmarkError(f"LinearRaccess output metadata mismatch: {mismatches}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--tier", required=True)
    parser.add_argument("--software-dir", type=Path, required=True)
    parser.add_argument("--linracc", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest_path = args.manifest.resolve()
        validated = validate_manifest(manifest_path, args.tier)
        paper_snapshot = git_snapshot(PAPER_ROOT)
        if paper_snapshot["status"]:
            raise BenchmarkError("paper repository must be clean before a benchmark run")
        lock = load_json_object(LOCK_PATH)
        software_dir = args.software_dir.resolve()
        software_snapshot = verify_software_checkout(software_dir, lock)
        executable = args.linracc.resolve()
        if not executable.is_file():
            raise BenchmarkError(f"LinearRaccess executable does not exist: {executable}")
        if not executable.stat().st_mode & 0o111:
            raise BenchmarkError(f"LinearRaccess executable is not executable: {executable}")

        build_info = verify_binary_build(executable, software_snapshot["commit"])
        output_dir = args.output_dir.resolve()
        if output_dir.exists():
            raise BenchmarkError(f"result bundle already exists: {output_dir}")
        output_dir.parent.mkdir(parents=True, exist_ok=True)

        manifest = validated["manifest"]
        source = manifest["source"]
        sequence = validated["sequence"]
        source_name = validated["source_name"]
        lengths = validated["lengths"]
        access_lens = validated["access_lens"]
        config = validated["configuration"]
        tier = validated["tier"]
        software_commit = software_snapshot["commit"]
        raw_rows: list[dict[str, Any]] = []
        summary_rows: list[dict[str, Any]] = []

        staging_raw = tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent)
        staging = Path(staging_raw)
        try:
            with tempfile.TemporaryDirectory(prefix="linearraccess_scaling_") as temp_raw:
                temporary = Path(temp_raw)
                for length in lengths:
                    fasta_path = temporary / f"prefix_{length}.fa"
                    fasta_path.write_text(
                        f">{source_name}_prefix_{length}\n{sequence[:length]}\n",
                        encoding="ascii",
                    )
                    result_path = temporary / f"prefix_{length}.txt"
                    command = command_for(
                        executable,
                        fasta_path,
                        result_path,
                        access_lens,
                        config,
                        software_commit,
                    )
                    for _ in range(tier["warmup_runs"]):
                        run_timed(command)

                    seconds: list[float] = []
                    rss_values: list[float] = []
                    output_hashes: list[str] = []
                    for repetition in range(1, tier["repetitions"] + 1):
                        elapsed, peak_rss_mb = run_timed(command)
                        metadata = parse_output_metadata(result_path)
                        check_metadata(metadata, access_lens, config, software_commit)
                        output_sha = sha256_file(result_path)
                        seconds.append(elapsed)
                        output_hashes.append(output_sha)
                        if peak_rss_mb is not None:
                            rss_values.append(peak_rss_mb)
                        raw_rows.append(
                            {
                                "tier": args.tier,
                                "source_id": source["id"],
                                "source_sequence_sha256": source["sequence_sha256"],
                                "prefix_sha256": sha256_text(sequence[:length]),
                                "length": length,
                                "repetition": repetition,
                                "seconds": elapsed,
                                "peak_rss_mb": peak_rss_mb,
                                "output_sha256": output_sha,
                                "configuration": config["id"],
                                "beam": config["beam"],
                                "c_multi": config["c_multi"],
                                "c_hairpin": config["c_hairpin"],
                                "fast_logsumexp": int(config["fast_logsumexp"]),
                            }
                        )
                        print(
                            f"{length}\trepetition={repetition}\t{elapsed:.6f}s",
                            flush=True,
                        )
                    if len(set(output_hashes)) != 1:
                        raise BenchmarkError(
                            f"fixed-input output was not deterministic for length {length}"
                        )
                    summary_rows.append(
                        {
                            "length": length,
                            "median_seconds": statistics.median(seconds),
                            "min_seconds": min(seconds),
                            "max_seconds": max(seconds),
                            "median_peak_rss_mb": (
                                statistics.median(rss_values) if rss_values else None
                            ),
                            "max_peak_rss_mb": max(rss_values) if rss_values else None,
                            "repetitions": tier["repetitions"],
                            "output_sha256": output_hashes[0],
                        }
                    )

            fits: dict[str, Any] = {
                "runtime_all_lengths": fit_log_log(summary_rows, "median_seconds"),
                "runtime_largest_three": fit_log_log(
                    summary_rows[-3:], "median_seconds"
                ),
            }
            rss_rows = [
                row for row in summary_rows if row["median_peak_rss_mb"] is not None
            ]
            if len(rss_rows) >= 2:
                fits["memory_all_lengths"] = fit_log_log(
                    rss_rows, "median_peak_rss_mb"
                )

            write_tsv(staging / "raw.tsv", raw_rows)
            write_tsv(staging / "summary.tsv", summary_rows)
            summary = {
                "schema_version": 1,
                "experiment_id": manifest["experiment_id"],
                "tier": args.tier,
                "measurements": summary_rows,
                "fits": fits,
            }
            (staging / "summary.json").write_text(
                json.dumps(summary, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            shutil.copyfile(manifest_path, staging / "manifest.json")
            command_template = command_for(
                executable,
                Path("{prefix_fasta}"),
                Path("{output}"),
                access_lens,
                config,
                software_commit,
            )
            provenance = {
                "schema_version": 1,
                "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                "experiment_id": manifest["experiment_id"],
                "tier": args.tier,
                "paper_repository": paper_snapshot,
                "software_repository": software_snapshot,
                "software_lock_sha256": sha256_file(LOCK_PATH),
                "manifest_path": str(manifest_path.relative_to(PAPER_ROOT)),
                "manifest_sha256": sha256_file(manifest_path),
                "source_path": str(validated["source_path"]),
                "source_file_sha256": source["file_sha256"],
                "source_sequence_sha256": source["sequence_sha256"],
                "linracc_path": str(executable),
                "linracc_sha256": sha256_file(executable),
                        "build_info": build_info,
                "command_template": command_template,
                "hardware": hardware_summary(),
                "argv": sys.argv,
                "timing": {
                    "clock": "time.perf_counter",
                    "peak_rss_source": "/usr/bin/time when available",
                },
            }
            (staging / "provenance.json").write_text(
                json.dumps(provenance, indent=2, sort_keys=True) + "\n",
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
