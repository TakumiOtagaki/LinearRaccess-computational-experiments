#!/usr/bin/env python3
"""Compare sampling-based window-accessibility estimates with direct
LinearRaccess posterior marginalization (pre-registered comparison)."""

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
from pathlib import Path
from typing import Any

import numpy as np

from benchmark_common import (
    BenchmarkError,
    git_snapshot,
    hardware_summary,
    load_json_object,
    read_single_fasta,
    sha256_file,
    sha256_text,
    verify_software_checkout,
    write_tsv,
)

PAPER_ROOT = Path(__file__).resolve().parents[2]
LOCK_PATH = PAPER_ROOT / "software.lock.json"
BINS = [1e-3, 1e-2, 1e-1, 0.5, 1.0]
BIN_LABELS = [
    "[0,1e-3)",
    "[1e-3,1e-2)",
    "[1e-2,1e-1)",
    "[1e-1,0.5)",
    "[0.5,1]",
]


def validate_manifest(path: Path, tier_name: str) -> dict[str, Any]:
    manifest = load_json_object(path)
    if manifest.get("schema_version") != 1:
        raise BenchmarkError("sampling-compare manifest must use schema_version 1")
    spec = manifest.get("sampling_software")
    if not isinstance(spec, dict) or not spec.get("commit"):
        raise BenchmarkError("sampling_software.commit required")
    window_lens = manifest.get("window_lens")
    if not isinstance(window_lens, list) or sorted(set(window_lens)) != window_lens:
        raise BenchmarkError("window_lens must be unique and sorted")

    sequences = manifest.get("sequences")
    if not isinstance(sequences, list):
        raise BenchmarkError("manifest requires a sequences list")
    for seq in sequences:
        src = Path(str(seq["path"]))
        if src.is_absolute() or not (path.parent / src).is_file():
            raise BenchmarkError(f"sequence FASTA does not exist: {src}")
        name, sequence = read_single_fasta(path.parent / src)
        if len(sequence) != seq.get("length"):
            raise BenchmarkError(f"length mismatch for {seq['id']}")
        if sha256_file(path.parent / src) != seq.get("file_sha256"):
            raise BenchmarkError(f"file SHA-256 mismatch for {seq['id']}")
        if sha256_text(sequence) != seq.get("sequence_sha256"):
            raise BenchmarkError(f"sequence SHA-256 mismatch for {seq['id']}")

    tiers = manifest.get("tiers")
    if not isinstance(tiers, dict) or tier_name not in tiers:
        raise BenchmarkError(f"unknown tier {tier_name!r}")
    tier = tiers[tier_name]
    for key in ("reference", "estimate"):
        if not isinstance(tier.get(key), dict):
            raise BenchmarkError(f"tier.{key} must be an object")
    for key in ("samples",):
        if tier["estimate"].get(key) is None and tier["estimate"].get("samples") is None:
            raise BenchmarkError("tier.estimate.samples required")
    return {"manifest": manifest, "sequences": sequences, "tier": tier}


def stream_samples(
    sampling_dir: Path,
    binary: str,
    beamsize: int,
    k: int,
    seq_file: Path,
    length: int,
):
    """Stream bracket-pair sample lines from a LinearSampling binary."""
    executable = sampling_dir / binary
    if not executable.is_file():
        raise BenchmarkError(f"LinearSampling binary missing: {executable}")
    command = [str(executable), str(beamsize), "0", "0", str(k), "", "0", "", ""]
    with seq_file.open("r", encoding="ascii") as handle:
        proc = subprocess.Popen(
            command, stdin=handle, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True,
        )
        assert proc.stdout is not None and proc.stderr is not None
        num_samples = 0
        for line in proc.stdout:
            line = line.strip()
            if len(line) != length or not set(line) <= set("()."):
                continue
            num_samples += 1
            yield line
        proc.stdout.close()
        stderr = proc.stderr.read()
        proc.stderr.close()
        returncode = proc.wait()
    if returncode != 0:
        raise BenchmarkError(
            f"LinearSampling failed ({returncode}): {command} "
            f"stderr={stderr[:1000]!r}"
        )
    if num_samples != k:
        raise BenchmarkError(
            f"LinearSampling returned {num_samples}/{k} samples for {binary} b={beamsize}"
        )


def window_frequencies(
    samples: "iterable[str]",
    length: int,
    window_lens: list[int],
    count_only: bool = False,
) -> dict[int, np.ndarray]:
    counts: dict[int, np.ndarray] = {
        w: np.zeros(length - w + 1, dtype=np.int64) for w in window_lens
    }
    for sample in samples:
        if len(sample) != length:
            raise BenchmarkError("sample length does not match sequence")
        paired = np.frombuffer(sample.encode("ascii"), dtype=np.uint8)
        paired = (paired != ord(".")).astype(np.int64)
        cs = np.concatenate(([0], np.cumsum(paired)))
        for w in window_lens:
            counts[w] += (cs[w:] - cs[:-w]) == 0
    return counts


def parse_linracc(path: Path, window_lens: list[int]) -> dict[int, np.ndarray]:
    probs: dict[int, np.ndarray] = {w: [] for w in window_lens}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or line.startswith(">") or not line.strip():
            continue
        head, _, tail = line.partition("\t")
        if not tail:
            continue
        start = int(head)
        for token in tail.split(";"):
            token = token.strip()
            if not token:
                continue
            w_str, _, p_str = token.partition(",")
            w = int(w_str)
            if w in probs:
                probs[w].append((start, float(p_str)))
    result: dict[int, np.ndarray] = {}
    for w, pairs in probs.items():
        pairs.sort()
        result[w] = np.array([p for _, p in pairs])
    return result


def run_linracc(
    linracc: Path,
    fasta_path: Path,
    output_path: Path,
    window_lens: list[int],
    beam: int,
    software_commit: str,
) -> dict[int, np.ndarray]:
    command = [
        str(linracc),
        f"-seqfile={fasta_path}",
        f"-outfile={output_path}",
        f"-access_len={','.join(map(str, window_lens))}",
        "-engine=lincapr",
        "-energy=turner2004",
        f"-beam={beam}",
        "-c_multi=30",
        "-c_hairpin=30",
        "-probabilities",
        "-metadata",
        f"-source-revision={software_commit}",
        "-no-normalize",
    ]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise BenchmarkError(
            f"LinearRaccess failed ({completed.returncode}): {' '.join(command)}\n"
            f"{completed.stderr[:1000]}"
        )
    return parse_linracc(output_path, window_lens)


def bin_index(p: float) -> int:
    for i, upper in enumerate(BINS):
        if p < upper:
            return i
    return len(BINS) - 1


def aggregate(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    aggregate_rows: list[dict[str, Any]] = []
    for row in rows:
        p_ref = row["p_ref"]
        if p_ref <= 0:
            continue
        for key, S in (("p100", 100), ("p1e3", 1000), ("p1e4", 10000)):
            estimates = row.get(key) or []
            if not estimates:
                continue
            rel_devs = [(float(e) - p_ref) / p_ref for e in estimates]
            rmse = math.sqrt(statistics.fmean(d * d for d in rel_devs))
            predicted = math.sqrt((1.0 - p_ref) / (S * p_ref))
            zero_frac = sum(1 for e in estimates if float(e) == 0.0) / len(estimates)
            aggregate_rows.append(
                {
                    "sequence": row["sequence"],
                    "w": row["w"],
                    "bin": BIN_LABELS[bin_index(p_ref)],
                    "p_ref": p_ref,
                    "S": S,
                    "median_rel_dev": statistics.median(rel_devs),
                    "rmse_rel_dev": rmse,
                    "predicted_rel_sd": predicted,
                    "zero_fraction": zero_frac,
                }
            )
    return aggregate_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--tier", required=True)
    parser.add_argument("--sampling-dir", type=Path, required=True)
    parser.add_argument("--software-dir", type=Path, required=True)
    parser.add_argument("--linracc", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest_path = args.manifest.resolve()
        validated = validate_manifest(manifest_path, args.tier)
        tiers = validated["manifest"]["tiers"]
        tier = validated["tier"]
        sequences = validated["sequences"]
        window_lens = validated["manifest"]["window_lens"]

        paper_snapshot = git_snapshot(PAPER_ROOT)
        if paper_snapshot["status"]:
            raise BenchmarkError("paper repository must be clean before a benchmark run")
        lock = load_json_object(LOCK_PATH)
        software_snapshot = verify_software_checkout(args.software_dir.resolve(), lock)
        software_commit = software_snapshot["commit"]

        sampling_dir = args.sampling_dir.resolve()
        sampling_snapshot = git_snapshot(sampling_dir)
        spec = validated["manifest"]["sampling_software"]
        if sampling_snapshot["commit"] != spec["commit"]:
            raise BenchmarkError(
                f"LinearSampling commit mismatch: expected {spec['commit']}, "
                f"got {sampling_snapshot['commit']}"
            )
        tracked_proc = subprocess.run(
            ["git", "diff", "--quiet", "HEAD"],
            cwd=sampling_dir,
            capture_output=True,
            check=False,
        )
        if tracked_proc.returncode != 0:
            raise BenchmarkError(
                f"LinearSampling checkout has tracked modifications: {sampling_dir}"
            )

        linracc_path = args.linracc.resolve()
        if not linracc_path.is_file() or not linracc_path.stat().st_mode & 0o111:
            raise BenchmarkError(f"LinearRaccess executable missing: {linracc_path}")

        if args.tier == "genome":
            genome_id = tier["sequence"]
            sequences = [s for s in sequences if s["id"] == genome_id]
            if not sequences:
                raise BenchmarkError(f"genome tier sequence not found: {genome_id}")
        elif tier.get("sequences"):
            allowed = set(tier["sequences"])
            sequences = [s for s in sequences if s["id"] in allowed]
            if not sequences:
                raise BenchmarkError("tier.sequences matched no manifest sequences")

        output_dir = args.output_dir.resolve()
        if output_dir.exists():
            raise BenchmarkError(f"result bundle already exists: {output_dir}")
        output_dir.parent.mkdir(parents=True, exist_ok=True)

        staging_raw = tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent)
        staging = Path(staging_raw)
        try:
            with tempfile.TemporaryDirectory(prefix="sampling_compare_") as temp_raw:
                temporary = Path(temp_raw)
                ref_mode = tier["reference"]["mode"]
                ref_samples = int(tier["reference"]["samples"])
                ref_beam = int(tier["reference"].get("beam", 0))
                bin_pick = {
                    "exact": spec["exact_binary"],
                    "linear": spec["linear_binary"],
                }[ref_mode]
                s_samples = [int(s) for s in tier["estimate"]["samples"]]
                s_beam = int(tier["estimate"]["beam"])
                replicates = int(tier["estimate"]["replicates"])

                raw_rows: list[dict[str, Any]] = []
                summary_rows: list[dict[str, Any]] = []
                for seq in sequences:
                    seq_path = (manifest_path.parent / seq["path"]).resolve()
                    _, sequence = read_single_fasta(seq_path)
                    length = len(sequence)
                    plain = temporary / f"{seq['id']}.seq"
                    plain.write_text(sequence + "\n", encoding="ascii")

                    ref_counts = window_frequencies(
                        stream_samples(
                            sampling_dir, bin_pick, ref_beam, ref_samples, plain, length
                        ),
                        length,
                        window_lens,
                    )
                    p_ref = {
                        w: ref_counts[w] / ref_samples for w in window_lens
                    }

                    estimates: dict[int, dict[int, list[list[float]]]] = {}
                    for S in s_samples:
                        estimates[S] = {
                            w: [[] for _ in range(length - w + 1)] for w in window_lens
                        }
                        for rep in range(1, replicates + 1):
                            freqs = window_frequencies(
                                stream_samples(
                                    sampling_dir,
                                    spec["linear_binary"],
                                    s_beam,
                                    S,
                                    plain,
                                    length,
                                ),
                                length,
                                window_lens,
                            )
                            for w in window_lens:
                                for start, value in enumerate(freqs[w].tolist()):
                                    estimates[S][w][start].append(float(value) / S)

                    p_star: dict[int, np.ndarray] | None = None
                    if args.tier == "genome":
                        linracc_out = temporary / f"{seq['id']}.linracc.txt"
                        beam = int(tier.get("linracc_beam", 200))
                        p_star = run_linracc(
                            linracc_path,
                            seq_path,
                            linracc_out,
                            window_lens,
                            beam,
                            software_commit,
                        )

                    for w in window_lens:
                        n_windows = length - w + 1
                        for start in range(n_windows):
                            row = {
                                "sequence": seq["id"],
                                "w": w,
                                "start": start,
                                "p_ref": float(p_ref[w][start]),
                                "p100": list(estimates.get(100, {}).get(w, [])[start]),
                                "p1e3": list(estimates.get(1000, {}).get(w, [])[start]),
                                "p1e4": list(estimates.get(10000, {}).get(w, [])[start]),
                                "p_star": (
                                    float(p_star[w][start]) if p_star is not None else None
                                ),
                            }
                            raw_rows.append(row)
                    summary_rows.extend(aggregate(raw_rows[-n_windows * len(window_lens):]))

                write_tsv(staging / "raw.tsv", raw_rows)
                if summary_rows:
                    write_tsv(staging / "aggregate.tsv", summary_rows)
                summary = {
                    "schema_version": 1,
                    "experiment_id": validated["manifest"]["experiment_id"],
                    "tier": args.tier,
                    "sequences": [s["id"] for s in sequences],
                    "window_lens": window_lens,
                    "reference": tier["reference"],
                    "estimate": tier["estimate"],
                    "rows": len(raw_rows),
                }
                (staging / "summary.json").write_text(
                    json.dumps(summary, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8",
                )
                shutil.copyfile(manifest_path, staging / "manifest.json")
                provenance = {
                    "schema_version": 1,
                    "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                    "experiment_id": validated["manifest"]["experiment_id"],
                    "tier": args.tier,
                    "paper_repository": paper_snapshot,
                    "software_repository": software_snapshot,
                    "software_lock_sha256": sha256_file(LOCK_PATH),
                    "sampling_repository": sampling_snapshot,
                    "manifest_path": str(manifest_path.relative_to(PAPER_ROOT)),
                    "manifest_sha256": sha256_file(manifest_path),
                    "linracc_path": str(linracc_path),
                    "linracc_sha256": sha256_file(linracc_path),
                    "linear_binary_sha256": sha256_file(sampling_dir / spec["linear_binary"]),
                    "exact_binary_sha256": sha256_file(sampling_dir / spec["exact_binary"]),
                    "hardware": hardware_summary(),
                    "argv": sys.argv,
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