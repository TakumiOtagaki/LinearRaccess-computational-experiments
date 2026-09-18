from __future__ import annotations

import math
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from benchmark_common import (  # noqa: E402
    BenchmarkError,
    fit_log_log,
    parse_output_metadata,
    read_single_fasta,
    sha256_text,
)
from run_scaling import check_metadata, command_for  # noqa: E402


class FastaTests(unittest.TestCase):
    def test_normalizes_dna_and_whitespace(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "input.fa"
            path.write_text(">id\nac gt\n", encoding="utf-8")
            name, sequence = read_single_fasta(path)
        self.assertEqual(name, "id")
        self.assertEqual(sequence, "ACGU")
        self.assertEqual(
            sha256_text(sequence),
            "a3c532aca6cd82dd0949dc252d6d1031da32eefbb7f2f42404923359b1a9f08a",
        )

    def test_rejects_multiple_records(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "input.fa"
            path.write_text(">one\nACGU\n>two\nACGU\n", encoding="utf-8")
            with self.assertRaisesRegex(BenchmarkError, "one record"):
                read_single_fasta(path)


class FitTests(unittest.TestCase):
    def test_recovers_quadratic_slope(self) -> None:
        rows = [
            {"length": 10, "value": 100},
            {"length": 20, "value": 400},
            {"length": 40, "value": 1600},
        ]
        fit = fit_log_log(rows, "value")
        self.assertTrue(math.isclose(fit["slope"], 2.0, abs_tol=1e-12))
        self.assertTrue(math.isclose(fit["r_squared"], 1.0, abs_tol=1e-12))


class MetadataTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = {
            "engine": "lincapr",
            "energy": "turner2004",
            "beam": 20,
            "c_multi": 30,
            "c_hairpin": 30,
            "fast_logsumexp": True,
            "normalize_profiles": False,
        }
        self.revision = "51ed3c3c5cf20e61cf075d116903baaf98cf4790"

    def test_command_contains_pinned_metadata(self) -> None:
        command = command_for(
            Path("LinRacc"),
            Path("input.fa"),
            Path("output.txt"),
            [1, 3],
            self.config,
            self.revision,
        )
        self.assertFalse(any(arg.startswith("-source-revision=") for arg in command))
        self.assertIn("-metadata", command)
        self.assertIn("-probabilities", command)
        self.assertIn("-no-normalize", command)

    def test_parses_and_checks_output_metadata(self) -> None:
        content = "\n".join(
            [
                "# linearraccess_metadata=1",
                f"# source_revision={self.revision}",
                "# engine=lincapr",
                "# effective_energy_model=turner2004",
                "# beam=20",
                "# c_multi=30",
                "# c_hairpin=30",
                "# logsumexp=fast",
                "# normalize_profiles=false",
                "# output_mode=probability",
                "# access_lens=1,3",
                ">id",
            ]
        )
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "output.txt"
            path.write_text(content + "\n", encoding="utf-8")
            metadata = parse_output_metadata(path)
        check_metadata(metadata, [1, 3], self.config, self.revision)

    def test_rejects_metadata_mismatch(self) -> None:
        with self.assertRaisesRegex(BenchmarkError, "metadata mismatch"):
            check_metadata({}, [1, 3], self.config, self.revision)


if __name__ == "__main__":
    unittest.main()
