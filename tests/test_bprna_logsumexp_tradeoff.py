"""Validate benchmark parsing and the agreed numerical tolerance."""
from pathlib import Path
import json
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from run_bprna_logsumexp_tradeoff import bounds, parse_output, summarize


class TradeoffTests(unittest.TestCase):
    def output(self):
        header = '\n'.join(f'# {k}={v}' for k, v in dict(
            beam=100, c_multi=30, c_hairpin=30, engine='lincapr',
            effective_energy_model='turner2004', output_mode='probability',
            access_lens='1,3,7,10', logsumexp='fast').items())
        record = ''.join(f"{i}\t" + ''.join(f'{k},0.5;' for k in (1, 3, 7, 10) if i+k <= 10) + '\n' for i in range(10))
        return header + '\n>a\n' + record + '>b\n' + record

    def parse(self, text, mode='fast'):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'output.txt'
            path.write_text(text)
            return parse_output(path, ['a', 'b'], 10, mode)

    def test_all_repeated_records_and_valid_terminal_windows(self):
        self.assertEqual(len(self.parse(self.output())), 10+8+4+1)

    def test_repetition_mismatch(self):
        with self.assertRaises(ValueError):
            self.parse(self.output().replace('1,0.5;', '1,0.4;', 1))

    def test_nonfinite_rejected(self):
        with self.assertRaises(ValueError):
            self.parse(self.output().replace('1,0.5;', '1,nan;', 1))

    def test_missing_window_rejected(self):
        with self.assertRaises(ValueError):
            self.parse(self.output().replace('10,0.5;', '', 1))

    def test_wrong_mode_rejected(self):
        with self.assertRaises(ValueError):
            self.parse(self.output(), 'exact')

    def test_tiny_negative_is_retained_but_tolerated(self):
        result = bounds([-1e-13, 0.5, 1.0001])
        self.assertEqual(result['min'], -1e-13)
        self.assertEqual(result['negative'], 1)
        self.assertEqual(result['outside_1e_12'], 1)

    def test_primary_aggregate_excludes_selected_anomalies(self):
        rows = []
        for name, cohort, speed, error in [('random', 'random', 1.2, .001),
                                          ('anomaly', 'sentinel', 4.0, .2)]:
            rows.append(dict(seq_id=name, cohort=cohort, stratum='50-99' if cohort == 'random' else 'sentinel',
                             length=50, batch_size=16, fast_median_sec=1.0, exact_median_sec=speed,
                             speedup_exact_over_fast=speed,
                             stats_by_lens={1: dict(abs_error={'max': error},
                                                   fast={'count': 50, 'outside_1e_12': 1},
                                                   exact={'count': 50, 'outside_1e_12': 0})}))
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)
            summarize(out, rows)
            result = json.loads((out / 'summary.json').read_text())['groups']['all_random']
            self.assertEqual(result['n'], 1)
            self.assertEqual(result['speedup']['median'], 1.2)
            self.assertEqual(result['max_abs_probability_difference'], .001)


if __name__ == '__main__':
    unittest.main()
