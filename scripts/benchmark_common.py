#!/usr/bin/env python3
"""Dependency-free helpers shared by LinearRaccess paper benchmarks."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


FULL_SHA = re.compile(r"^[0-9a-f]{40}$")


class BenchmarkError(RuntimeError):
    """A benchmark precondition or execution invariant failed."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("ascii")).hexdigest()


def load_json_object(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise BenchmarkError(f"cannot read valid JSON from {path}: {error}") from error
    if not isinstance(value, dict):
        raise BenchmarkError(f"JSON root must be an object: {path}")
    return value


def read_single_fasta(path: Path) -> tuple[str, str]:
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
                raise BenchmarkError(f"benchmark FASTA must contain one record: {path}")
            name = line[1:].strip()
            if not name:
                raise BenchmarkError(f"FASTA header must not be empty: {path}")
        elif name is None:
            raise BenchmarkError(f"sequence occurs before the FASTA header: {path}")
        else:
            chunks.append("".join(line.split()).upper().replace("T", "U"))

    if name is None or not chunks:
        raise BenchmarkError(f"FASTA contains no sequence: {path}")
    sequence = "".join(chunks)
    invalid = next((index for index, base in enumerate(sequence) if base not in "ACGU"), None)
    if invalid is not None:
        raise BenchmarkError(
            f"unsupported nucleotide at 1-based position {invalid + 1} in {path}"
        )
    return name, sequence


def git_output(args: list[str], directory: Path) -> str:
    completed = subprocess.run(
        ["git", *args],
        cwd=directory,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise BenchmarkError(f"git {' '.join(args)} failed in {directory}: {detail}")
    return completed.stdout.strip()


def git_snapshot(directory: Path) -> dict[str, str]:
    return {
        "commit": git_output(["rev-parse", "HEAD"], directory),
        "status": git_output(["status", "--porcelain=v1"], directory),
        "remote": git_output(["remote", "get-url", "origin"], directory),
    }


def normalize_github_remote(remote: str) -> str:
    value = remote.removesuffix(".git")
    prefix = "git@github.com:"
    if value.startswith(prefix):
        return f"https://github.com/{value[len(prefix):]}"
    return value


def verify_software_checkout(
    directory: Path, lock: dict[str, Any]
) -> dict[str, str]:
    if not directory.is_dir():
        raise BenchmarkError(f"software checkout does not exist: {directory}")
    software = lock.get("software")
    if not isinstance(software, dict):
        raise BenchmarkError("software.lock.json requires a software object")
    expected_commit = str(software.get("commit", ""))
    if not FULL_SHA.fullmatch(expected_commit):
        raise BenchmarkError("software lock does not contain a full lowercase Git SHA")

    snapshot = git_snapshot(directory)
    if snapshot["commit"] != expected_commit:
        raise BenchmarkError(
            f"software commit mismatch: expected {expected_commit}, "
            f"got {snapshot['commit']}"
        )
    if snapshot["status"]:
        raise BenchmarkError(f"software checkout is dirty: {directory}")
    expected_remote = str(software.get("repository", ""))
    actual_remote = normalize_github_remote(snapshot["remote"])
    if actual_remote != expected_remote:
        raise BenchmarkError(
            f"software remote mismatch: expected {expected_remote}, got {actual_remote}"
        )
    return snapshot


def verify_binary_build(executable: Path, expected_revision: str) -> dict[str, str]:
    """Check the revision compiled into the executable, without a CLI override."""
    try:
        result = subprocess.run([str(executable), "--build-info"], check=True,
                                capture_output=True, text=True, timeout=30)
        info = json.loads(result.stdout)
    except (OSError, subprocess.SubprocessError, ValueError) as error:
        raise BenchmarkError(f"cannot read build information from {executable}: {error}") from error
    if (not isinstance(info, dict)
            or not FULL_SHA.fullmatch(str(info.get("source_revision", "")))
            or info["source_revision"] != expected_revision
            or not isinstance(info.get("software_version"), str)):
        raise BenchmarkError(f"binary build does not match clean revision {expected_revision}: {info}")
    return info


def run_timed(command: list[str]) -> tuple[float, float | None]:
    time_executable = Path("/usr/bin/time")
    wrapped = command
    parser: str | None = None
    if time_executable.exists():
        if sys.platform == "darwin":
            wrapped = [str(time_executable), "-l", *command]
            parser = "darwin"
        elif sys.platform.startswith("linux"):
            wrapped = [str(time_executable), "-v", *command]
            parser = "linux"

    started = time.perf_counter()
    completed = subprocess.run(wrapped, text=True, capture_output=True, check=False)
    elapsed = time.perf_counter() - started
    if (
        completed.returncode != 0
        and parser == "darwin"
        and "sysctl kern.clockrate: Operation not permitted" in completed.stderr
    ):
        parser = None
        started = time.perf_counter()
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
        elapsed = time.perf_counter() - started
    if completed.returncode != 0:
        raise BenchmarkError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )

    peak_rss_mb: float | None = None
    if parser == "darwin":
        match = re.search(
            r"(?m)^\s*(\d+)\s+maximum resident set size\s*$", completed.stderr
        )
        if match:
            peak_rss_mb = int(match.group(1)) / (1024 * 1024)
    elif parser == "linux":
        match = re.search(
            r"Maximum resident set size \(kbytes\):\s*(\d+)", completed.stderr
        )
        if match:
            peak_rss_mb = int(match.group(1)) / 1024
    return elapsed, peak_rss_mb


def fit_log_log(rows: list[dict[str, Any]], value_key: str) -> dict[str, float]:
    if len(rows) < 2:
        raise BenchmarkError("a log-log fit requires at least two observations")
    points: list[tuple[float, float]] = []
    for row in rows:
        x = float(row["length"])
        y = float(row[value_key])
        if x <= 0.0 or y <= 0.0 or not (math.isfinite(x) and math.isfinite(y)):
            raise BenchmarkError(f"log-log fit requires positive finite values: {row}")
        points.append((math.log(x), math.log(y)))

    x_mean = statistics.fmean(x for x, _ in points)
    y_mean = statistics.fmean(y for _, y in points)
    denominator = sum((x - x_mean) ** 2 for x, _ in points)
    if denominator == 0.0:
        raise BenchmarkError("log-log fit requires distinct sequence lengths")
    slope = sum((x - x_mean) * (y - y_mean) for x, y in points) / denominator
    intercept = y_mean - slope * x_mean
    fitted = [intercept + slope * x for x, _ in points]
    residual_ss = sum((y - prediction) ** 2 for (_, y), prediction in zip(points, fitted))
    total_ss = sum((y - y_mean) ** 2 for _, y in points)
    return {
        "slope": slope,
        "intercept": intercept,
        "r_squared": 1.0 - residual_ss / total_ss if total_ss else 1.0,
    }


def parse_output_metadata(path: Path) -> dict[str, str]:
    metadata: dict[str, str] = {}
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise BenchmarkError(f"cannot read LinearRaccess output {path}: {error}") from error
    for line in lines:
        if not line.startswith("# "):
            continue
        key, separator, value = line[2:].partition("=")
        if not separator or not key:
            continue
        if key in metadata:
            raise BenchmarkError(f"duplicate output metadata key {key!r} in {path}")
        metadata[key] = value
    return metadata


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise BenchmarkError(f"refusing to write an empty table: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def hardware_summary() -> dict[str, str]:
    cpu = platform.processor()
    if not cpu and sys.platform.startswith("linux"):
        try:
            for line in Path("/proc/cpuinfo").read_text(encoding="utf-8").splitlines():
                if line.lower().startswith("model name"):
                    cpu = line.partition(":")[2].strip()
                    break
        except OSError:
            pass
    if not cpu and sys.platform == "darwin":
        completed = subprocess.run(
            ["sysctl", "-n", "machdep.cpu.brand_string"],
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode == 0:
            cpu = completed.stdout.strip()
    return {
        "hostname": platform.node(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": cpu or "unknown",
        "python": sys.version,
    }
