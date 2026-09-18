import json
import sys
import subprocess
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from benchmark_common import BenchmarkError, verify_binary_build


class BuildProvenanceTests(unittest.TestCase):
    def test_clean_matching_build(self):
        info = {"source_revision": "a" * 40, "software_version": "0.1.0"}
        with patch("benchmark_common.subprocess.run", return_value=subprocess.CompletedProcess([], 0, json.dumps(info))):
            self.assertEqual(verify_binary_build(Path("/example/LinRacc"), "a" * 40), info)

    def test_rejects_unverifiable_or_mismatched_build(self):
        for revision in ["unknown", "a" * 40 + "-dirty", "b" * 40, None]:
            with self.subTest(revision=revision), patch("benchmark_common.subprocess.run", return_value=subprocess.CompletedProcess([], 0, json.dumps({"source_revision": revision, "software_version": "0.1.0"}))):
                with self.assertRaises(BenchmarkError):
                    verify_binary_build(Path("/example/LinRacc"), "a" * 40)

    def test_rejects_old_cli_and_malformed_json(self):
        with patch("benchmark_common.subprocess.run", side_effect=subprocess.CalledProcessError(1, [])):
            with self.assertRaises(BenchmarkError):
                verify_binary_build(Path("/example/LinRacc"), "a" * 40)
        with patch("benchmark_common.subprocess.run", return_value=subprocess.CompletedProcess([], 0, "not JSON")):
            with self.assertRaises(BenchmarkError):
                verify_binary_build(Path("/example/LinRacc"), "a" * 40)
