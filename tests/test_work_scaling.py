from __future__ import annotations

import sys
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from benchmark_common import BenchmarkError  # noqa: E402
from run_work_scaling import TABLES, parse_debug_work  # noqa: E402


def debug_log(mb_counts: tuple[int, int] = (2, 3)) -> str:
    lines: list[str] = []
    for j in range(2):
        for table in TABLES:
            lines.append(
                f"debug_beam: j={j} table={table} before=3 after=2 threshold=0"
            )
        lines.append(
            f"debug_beam_updates: j={j} M2={j + 1} MB={mb_counts[j]} "
            f"M1={j + 3} M={j + 4}"
        )
    return "\n".join(lines)


class WorkCounterTests(unittest.TestCase):
    def test_parses_counts_and_bound(self) -> None:
        row = parse_debug_work(debug_log(), length=2, beam=2)
        self.assertEqual(row["inside_m2_updates"], 3)
        self.assertEqual(row["inside_mb_bifurcation_updates"], 5)
        self.assertEqual(row["mb_bifurcation_n_b2_bound"], 8)
        self.assertEqual(row["mb_bifurcation_bound_fraction"], 5 / 8)
        self.assertEqual(row["m2_max_after_prune"], 2)

    def test_rejects_missing_position(self) -> None:
        incomplete = "\n".join(debug_log().splitlines()[:-1])
        with self.assertRaises(BenchmarkError):
            parse_debug_work(incomplete, length=2, beam=2)

    def test_rejects_count_above_bound(self) -> None:
        with self.assertRaisesRegex(BenchmarkError, "N\\*b\\^2 bound"):
            parse_debug_work(debug_log((5, 5)), length=2, beam=2)


if __name__ == "__main__":
    unittest.main()
