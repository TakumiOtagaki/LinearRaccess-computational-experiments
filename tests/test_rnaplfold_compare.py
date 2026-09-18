import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import rnaplfold_compare as comparison


class RNAplfoldComparisonTests(unittest.TestCase):
    def test_end_anchored_rows_become_start_anchored_windows(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "example_lunp"
            path.write_text("# unpaired\n1\t0.9\tNA\tNA\n3\t0.7\t0.6\t0.5\n")
            parsed = comparison.parse_lunp(path)
        self.assertEqual(parsed[0], {1: 0.9, 3: 0.5})
        self.assertEqual(parsed[1][2], 0.6)
        self.assertEqual(parsed[2][1], 0.7)

    def maps(self):
        exact, approx = {}, {}
        for length in (1, 3, 7, 10):
            for start in range(11 - length):
                exact.setdefault(start, {})[length] = 0.0
                approx.setdefault(start, {})[length] = 0.75
        return exact, approx

    def test_all_requested_cells_including_last_start_are_compared(self):
        exact, approx = self.maps()
        self.assertEqual(comparison.all_errors(exact, approx, 10), [0.25] * 23)

    def test_missing_boundary_cell_fails_instead_of_biasing_median(self):
        exact, approx = self.maps()
        del approx[0][10]
        with self.assertRaisesRegex(ValueError, "Missing cell"):
            comparison.all_errors(exact, approx, 10)

    def test_nonfinite_probability_fails(self):
        exact, approx = self.maps()
        approx[0][10] = float("nan")
        with self.assertRaisesRegex(ValueError, "Non-finite"):
            comparison.all_errors(exact, approx, 10)


if __name__ == "__main__":
    unittest.main()
