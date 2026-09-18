from __future__ import annotations

import argparse
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import run_raccess_compare as runner
from benchmark_common import BenchmarkError


class ComparisonBundleTests(unittest.TestCase):
    def arguments(self, root: Path) -> argparse.Namespace:
        source = root / "input.fa"
        source.write_text(">example\nAC\n", encoding="utf-8")
        return argparse.Namespace(
            bundle=root / "bundle", input=source, max_span=0,
            access_len="1", raccess_binary=source, linracc_binary=source,
            energy_thr=1000.0, linracc_engine="raccess", energy="turner1999",
            beam=0, source_revision=None, raccess_repo=root,
            raccess_commit="test", raccess_patch_note="",
        )

    def run_comparison(self, root: Path, observed: str) -> tuple[int, dict]:
        args = self.arguments(root)

        def fake_run(command):
            output = Path(next(x.split("=", 1)[1] for x in command if x.startswith("-outfile=")))
            content = observed if "-probabilities" in command else ">example\n0\t1,0;\n1\t1,0;\n"
            output.write_text(content, encoding="utf-8")
            return 0.01, 1.0

        with patch.object(runner, "parse_args", return_value=args), \
             patch.object(runner, "run_timed", side_effect=fake_run), \
             patch.object(runner, "verify_binary_build", return_value={"source_revision": "a" * 40, "software_version": "test"}), \
             patch.object(runner, "hardware_summary", return_value={}), \
             patch.object(runner, "git_sha", return_value="test"), \
             contextlib.redirect_stdout(io.StringIO()), \
             contextlib.redirect_stderr(io.StringIO()):
            status = runner.main()
        summary = json.loads((args.bundle / "summary.json").read_text())
        self.assertEqual(summary, json.loads((args.bundle / "analysis/summary.json").read_text()))
        self.assertTrue((args.bundle / "provenance.json").is_file())
        return status, summary

    def test_complete_comparison_writes_repository_compatible_summary(self):
        with tempfile.TemporaryDirectory() as raw:
            status, summary = self.run_comparison(Path(raw), ">example\n0\t1,1;\n1\t1,1;\n")
        self.assertEqual(status, 0)
        self.assertTrue(summary["comparison_complete"])
        self.assertEqual(summary["compared_cells"], 2)

    def test_missing_output_fails_but_keeps_diagnostics(self):
        with tempfile.TemporaryDirectory() as raw:
            status, summary = self.run_comparison(Path(raw), ">example\n0\t1,1;\n")
        self.assertEqual(status, 1)
        self.assertFalse(summary["comparison_complete"])
        self.assertEqual(summary["per_sequence"][0]["missing_linracc_cells"], 1)

    def test_nonfinite_energy_fails_but_keeps_diagnostics(self):
        with tempfile.TemporaryDirectory() as raw:
            status, summary = self.run_comparison(Path(raw), ">example\n0\t1,1;\n1\t1,0;\n")
        self.assertEqual(status, 1)
        self.assertEqual(summary["per_sequence"][0]["nonfinite_cells"], 1)

    def test_existing_bundle_is_not_modified(self):
        with tempfile.TemporaryDirectory() as raw:
            args = self.arguments(Path(raw))
            args.bundle.mkdir()
            sentinel = args.bundle / "summary.json"
            sentinel.write_text("original", encoding="utf-8")
            with patch.object(runner, "parse_args", return_value=args), \
                 patch.object(runner, "run_timed") as run:
                with self.assertRaisesRegex(BenchmarkError, "already exists"):
                    runner.main()
            run.assert_not_called()
            self.assertEqual(sentinel.read_text(), "original")


if __name__ == "__main__":
    unittest.main()
