#!/usr/bin/env python3
"""Count deterministic LinearRaccess inside-pass work across sequence prefixes."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import re
import shutil
import statistics
import subprocess
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
    sha256_file,
    sha256_text,
    verify_software_checkout,
    write_tsv,
)
from run_scaling import LOCK_PATH, PAPER_ROOT, validate_manifest


TABLES = ("S", "SE", "M", "MB", "M1", "M2")
UPDATE_RE = re.compile(
    r"^debug_beam_updates: j=(\d+) M2=(\d+) MB=(\d+) M1=(\d+) M=(\d+)$"
)
PRUNE_RE = re.compile(
    r"^debug_beam: j=(\d+) table=(S|SE|M|MB|M1|M2) "
    r"before=(\d+) after=(\d+) threshold="
)


def parse_debug_work(stderr: str, length: int, beam: int) -> dict[str, Any]:
    updates: dict[int, tuple[int, int, int, int]] = {}
    tables: dict[str, dict[int, tuple[int, int]]] = {table: {} for table in TABLES}
    for line in stderr.splitlines():
        update_match = UPDATE_RE.match(line)
        if update_match:
            j = int(update_match.group(1))
            if j in updates:
                raise BenchmarkError(f"duplicate update record for j={j}")
            updates[j] = tuple(int(update_match.group(index)) for index in range(2, 6))
            continue
        prune_match = PRUNE_RE.match(line)
        if prune_match:
            j = int(prune_match.group(1))
            table = prune_match.group(2)
            if j in tables[table]:
                raise BenchmarkError(f"duplicate {table} prune record for j={j}")
            tables[table][j] = (
                int(prune_match.group(3)),
                int(prune_match.group(4)),
            )

    expected_positions = set(range(length))
    if set(updates) != expected_positions:
        raise BenchmarkError(
            f"update records do not cover 0..{length - 1}: found {len(updates)}"
        )
    for table in TABLES:
        if set(tables[table]) != expected_positions:
            raise BenchmarkError(
                f"{table} prune records do not cover 0..{length - 1}: "
                f"found {len(tables[table])}"
            )
        if max(after for _, after in tables[table].values()) > beam:
            raise BenchmarkError(f"{table} retained more than beam={beam} states")

    m2_updates = sum(value[0] for value in updates.values())
    mb_bifurcation_updates = sum(value[1] for value in updates.values())
    m1_updates = sum(value[2] for value in updates.values())
    m_updates = sum(value[3] for value in updates.values())
    theoretical_bound = length * beam * beam
    if mb_bifurcation_updates > theoretical_bound:
        raise BenchmarkError(
            "observed MB bifurcation candidates exceed the fixed-beam N*b^2 bound"
        )

    row: dict[str, Any] = {
        "length": length,
        "inside_m2_updates": m2_updates,
        "inside_mb_bifurcation_updates": mb_bifurcation_updates,
        "inside_m1_updates": m1_updates,
        "inside_m_updates": m_updates,
        "inside_tracked_updates": (
            m2_updates + mb_bifurcation_updates + m1_updates + m_updates
        ),
        "mb_bifurcation_updates_per_nt": mb_bifurcation_updates / length,
        "mb_bifurcation_n_b2_bound": theoretical_bound,
        "mb_bifurcation_bound_fraction": (
            mb_bifurcation_updates / theoretical_bound
        ),
    }
    for table in TABLES:
        before = [tables[table][j][0] for j in range(length)]
        after = [tables[table][j][1] for j in range(length)]
        prefix = table.lower()
        row[f"{prefix}_mean_before_prune"] = statistics.fmean(before)
        row[f"{prefix}_mean_after_prune"] = statistics.fmean(after)
        row[f"{prefix}_max_before_prune"] = max(before)
        row[f"{prefix}_max_after_prune"] = max(after)
    return row


def command_for_work(
    executable: Path,
    fasta_path: Path,
    output_path: Path,
    access_len: int,
    config: dict[str, Any],
    length: int,
) -> list[str]:
    command = [
        str(executable),
        f"-seqfile={fasta_path}",
        f"-outfile={output_path}",
        f"-access_len={access_len}",
        f"-engine={config['engine']}",
        f"-beam={config['beam']}",
        f"-c_multi={config['c_multi']}",
        f"-c_hairpin={config['c_hairpin']}",
        "-debug-beam",
        f"0,{length - 1}",
        "-debug-beam-interval",
        "1",
        "-debug-skip-access",
    ]
    if config["engine"] == "lincapr":
        command.append(f"-energy={config['energy']}")
    if not config["normalize_profiles"]:
        command.append("-no-normalize")
    if not config["fast_logsumexp"]:
        command.append("-no-fast-logsumexp")
    return command


def pairwise_slopes(rows: list[dict[str, Any]], value_key: str) -> list[dict[str, Any]]:
    slopes: list[dict[str, Any]] = []
    for left, right in zip(rows, rows[1:]):
        left_value = float(left[value_key])
        right_value = float(right[value_key])
        if left_value <= 0.0 or right_value <= 0.0:
            raise BenchmarkError("pairwise log slope requires positive work counts")
        slopes.append(
            {
                "length_from": int(left["length"]),
                "length_to": int(right["length"]),
                "slope": math.log(right_value / left_value)
                / math.log(float(right["length"]) / float(left["length"])),
            }
        )
    return slopes


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
        software_snapshot = verify_software_checkout(args.software_dir.resolve(), lock)
        executable = args.linracc.resolve()
        if not executable.is_file():
            raise BenchmarkError(f"LinearRaccess executable does not exist: {executable}")
        if not executable.stat().st_mode & 0o111:
            raise BenchmarkError(f"LinearRaccess executable is not executable: {executable}")

        output_dir = args.output_dir.resolve()
        if output_dir.exists():
            raise BenchmarkError(f"result bundle already exists: {output_dir}")
        output_dir.parent.mkdir(parents=True, exist_ok=True)

        manifest = validated["manifest"]
        source = manifest["source"]
        sequence = validated["sequence"]
        source_name = validated["source_name"]
        lengths = validated["lengths"]
        access_len = validated["access_lens"][0]
        config = validated["configuration"]
        beam = int(config["beam"])
        rows: list[dict[str, Any]] = []

        staging_raw = tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent)
        staging = Path(staging_raw)
        try:
            with tempfile.TemporaryDirectory(prefix="linearraccess_work_") as temp_raw:
                temporary = Path(temp_raw)
                for length in lengths:
                    fasta_path = temporary / f"prefix_{length}.fa"
                    fasta_path.write_text(
                        f">{source_name}_prefix_{length}\n{sequence[:length]}\n",
                        encoding="ascii",
                    )
                    result_path = temporary / f"prefix_{length}.txt"
                    command = command_for_work(
                        executable,
                        fasta_path,
                        result_path,
                        access_len,
                        config,
                        length,
                    )
                    completed = subprocess.run(
                        command,
                        text=True,
                        capture_output=True,
                        check=False,
                    )
                    if completed.returncode != 0:
                        raise BenchmarkError(
                            f"command failed ({completed.returncode}): {' '.join(command)}\n"
                            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
                        )
                    row = parse_debug_work(completed.stderr, length, beam)
                    row = {
                        "tier": args.tier,
                        "source_id": source["id"],
                        "source_sequence_sha256": source["sequence_sha256"],
                        "prefix_sha256": sha256_text(sequence[:length]),
                        "configuration": config["id"],
                        "beam": beam,
                        "c_multi": config["c_multi"],
                        "c_hairpin": config["c_hairpin"],
                        **row,
                    }
                    rows.append(row)
                    print(
                        f"{length}\tmb_bifurcation_updates="
                        f"{row['inside_mb_bifurcation_updates']}\tbound_fraction="
                        f"{row['mb_bifurcation_bound_fraction']:.6f}",
                        flush=True,
                    )

            value_key = "inside_mb_bifurcation_updates"
            fits = {
                "mb_bifurcation_all_lengths": fit_log_log(rows, value_key),
                "mb_bifurcation_largest_three": fit_log_log(rows[-3:], value_key),
                "inside_tracked_all_lengths": fit_log_log(
                    rows, "inside_tracked_updates"
                ),
                "inside_tracked_largest_three": fit_log_log(
                    rows[-3:], "inside_tracked_updates"
                ),
                "mb_bifurcation_pairwise_slopes": pairwise_slopes(rows, value_key),
            }
            write_tsv(staging / "summary.tsv", rows)
            summary = {
                "schema_version": 1,
                "experiment_id": manifest["experiment_id"],
                "measurement": "inside-pass deterministic work counters",
                "tier": args.tier,
                "measurements": rows,
                "fits": fits,
            }
            (staging / "summary.json").write_text(
                json.dumps(summary, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            shutil.copyfile(manifest_path, staging / "manifest.json")
            command_template = command_for_work(
                executable,
                Path("{prefix_fasta}"),
                Path("{output}"),
                access_len,
                config,
                int(lengths[-1]),
            )
            command_template[command_template.index(f"0,{lengths[-1] - 1}")] = "0,{N-1}"
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
                "command_template": command_template,
                "hardware": hardware_summary(),
                "argv": sys.argv,
                "counter_source": "LinRacc -debug-beam with interval 1",
                "counter_scope": "inside-pass transitions and post-prune occupancy",
                "bound_audited": "MB update candidates <= N * beam^2",
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
