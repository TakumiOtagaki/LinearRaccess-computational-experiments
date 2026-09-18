from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from benchmark_common import BenchmarkError  # noqa: E402
from run_algorithm_validation import (  # noqa: E402
    average_ranks,
    invariant_counts,
    parse_probabilities,
    percentile,
    summarize_errors,
)


class ProbabilityParserTests(unittest.TestCase):
    def test_parses_complete_probability_grid(self) -> None:
        content = "\n".join(
            [
                "# output_mode=probability",
                ">example",
                "0\t1,0.9;2,0.7;",
                "1\t1,0.8;2,0.6;",
                "2\t1,0.7;",
            ]
        )
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "output.txt"
            path.write_text(content + "\n", encoding="utf-8")
            values = parse_probabilities(path, 3, [1, 2])
        self.assertEqual(values[(1, 2)], 0.6)
        self.assertEqual(len(values), 5)

    def test_rejects_incomplete_probability_grid(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "output.txt"
            path.write_text(">example\n0\t1,0.9;\n", encoding="utf-8")
            with self.assertRaisesRegex(BenchmarkError, "grid mismatch"):
                parse_probabilities(path, 2, [1])


class MetricTests(unittest.TestCase):
    def test_percentile_interpolates(self) -> None:
        self.assertEqual(percentile([0.0, 10.0], 50.0), 5.0)

    def test_average_ranks_handles_ties(self) -> None:
        self.assertEqual(average_ranks([10.0, 5.0, 10.0]), [2.5, 1.0, 2.5])

    def test_error_summary_and_invariants(self) -> None:
        reference = {(0, 1): 0.9, (1, 1): 0.8, (0, 2): 0.7}
        observed = {(0, 1): 0.8, (1, 1): 0.7, (0, 2): 0.6}
        summary = summarize_errors(observed, reference, 1e-8)
        self.assertAlmostEqual(summary["maximum_absolute_error"], 0.1)
        self.assertAlmostEqual(summary["spearman_rho"], 1.0)
        counts = invariant_counts(observed, 2, [1, 2], 1e-10)
        self.assertEqual(counts["probability_bound_violations"], 0)
        self.assertEqual(counts["nested_window_monotonicity_violations"], 0)


if __name__ == "__main__":
    unittest.main()
