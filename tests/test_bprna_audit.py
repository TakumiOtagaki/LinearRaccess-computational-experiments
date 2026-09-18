"""Regression checks for incomplete and corrupt probability outputs."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location(
    'audit_bprna', Path(__file__).resolve().parents[1] / 'scripts/audit_bprna.py')
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


class ProbabilityAuditTests(unittest.TestCase):
    def check_output(self, transform=lambda text: text):
        lines = [f'# {k}={v}' for k, v in audit.EXPECTED.items()]
        lines.append('>test')
        for pos in range(12):
            payload = ''.join(f'{k},0.5;' for k in audit.LENS if pos + k <= 12)
            lines.append(f'{pos}\t{payload}')
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / 'output.txt'
            p.write_text(transform('\n'.join(lines) + '\n'))
            return audit.audit_output(p, 'test', 12)

    def test_valid_terminal_windows(self):
        errors, _, stats = self.check_output()
        self.assertFalse(errors)
        self.assertEqual([stats[k]['count'] for k in audit.LENS], [12, 10, 6, 3])

    def test_nonfinite_and_out_of_range(self):
        for value, issue in [('nan', 'nonfinite_probability'),
                             ('inf', 'nonfinite_probability'),
                             ('-0.1', 'probability_out_of_range'),
                             ('1.01', 'probability_out_of_range')]:
            with self.subTest(value=value):
                errors, _, _ = self.check_output(lambda s: s.replace('1,0.5;', f'1,{value};', 1))
                self.assertEqual(errors[issue], 1)

    def test_missing_window(self):
        errors, _, _ = self.check_output(lambda s: s.replace('10,0.5;', '', 1))
        self.assertEqual(errors['window_coverage'], 1)

    def test_tolerance_diagnostics_do_not_hide_strict_failures(self):
        errors, _, stats = self.check_output(lambda s: s.replace('1,0.5;', '1,-1e-13;', 1))
        self.assertEqual(errors['probability_out_of_range'], 1)
        self.assertEqual(stats[1]['negative'], 1)
        self.assertEqual(stats[1]['outside_tolerance_1e_12'], 0)
        errors, _, stats = self.check_output(lambda s: s.replace('1,0.5;', '1,1.0008;', 1))
        self.assertEqual(errors['probability_out_of_range'], 1)
        self.assertEqual(stats[1]['above_one'], 1)
        self.assertEqual(stats[1]['outside_tolerance_1e_6'], 1)

    def test_truncated_output(self):
        errors, _, _ = self.check_output(lambda s: s.rsplit('11\t', 1)[0])
        self.assertEqual(errors['position_count'], 1)

    def test_wrong_identity_and_configuration(self):
        errors, _, _ = self.check_output(lambda s: s.replace('>test', '>other').replace('beam=100', 'beam=0'))
        self.assertEqual(errors['sequence_header'], 1)
        self.assertEqual(errors['metadata_beam'], 1)


if __name__ == '__main__':
    unittest.main()
