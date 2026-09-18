"""End-to-end checks for production task failures and output provenance."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / 'benchmark/scripts/run_bprna_chunk.py'
FAKE = r'''#!/usr/bin/env python3
import json, sys, time
from pathlib import Path
revision = 'a' * 40
if sys.argv[1:] == ['--build-info']:
    print(json.dumps(dict(software_version='0.1.0', source_revision=revision)))
    raise SystemExit(0)
options = dict(arg.split('=', 1) for arg in sys.argv[1:] if '=' in arg)
lines = Path(options['-seqfile']).read_text().splitlines()
name, sequence = lines[0][1:], ''.join(lines[1:])
if name == 'exit_fail': raise SystemExit(7)
if name == 'missing': raise SystemExit(0)
if name == 'timeout': time.sleep(5)
metadata = dict(software_version='0.1.0', source_revision=revision,
                engine='lincapr', effective_energy_model='turner2004',
                beam='100', c_multi='30', c_hairpin='30', access_lens='1,3,7,10',
                logsumexp='exact' if '-no-fast-logsumexp' in sys.argv else 'fast',
                output_mode='probability', normalize_profiles='true')
if name == 'wrong_revision': metadata['source_revision'] = 'b' * 40
with Path(options['-outfile']).open('w') as out:
    for k, v in metadata.items(): out.write(f'# {k}={v}\n')
    out.write('>' + name + '\n')
    for i in range(len(sequence) - (name == 'truncated')):
        values = []
        for k in (1,3,7,10):
            if i+k <= len(sequence):
                v = 'nan' if name == 'nan' and i == 0 else '0.5'
                if name == 'range' and i == 0: v = '1.0008'
                if name == 'range' and i == 1: v = '-1e-13'
                values.append(f'{k},{v};')
        out.write(str(i) + '\t' + ''.join(values) + '\n')
'''


class ChunkTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.binary = self.root / 'LinRacc'
        self.binary.write_text(FAKE)
        self.binary.chmod(0o755)
        self.chunk = self.root / 'chunk_01.fa'
        self.output = self.root / 'results'

    def launch(self, names, timeout=10, wrapper=False):
        self.chunk.write_text(''.join(f'>{name}\n' + 'A'*50 + '\n' for name in names))
        if wrapper:
            env = dict(os.environ, BPRNA_REPO_ROOT=str(ROOT), BPRNA_CHUNK_DIR=str(self.root),
                       BPRNA_OUTDIR=str(self.output), BPRNA_LINRACC=str(self.binary),
                       SGE_TASK_ID='1', JOB_ID='test-job', BPRNA_PARALLEL='2')
            return subprocess.run(['bash', str(ROOT / 'benchmark/scripts/uge/bprna_chunk_job.sh')],
                                  env=env, capture_output=True, text=True)
        return subprocess.run(['python3', str(RUNNER), '--chunk', str(self.chunk),
                               '--outdir', str(self.output), '--linracc', str(self.binary),
                               '--task-id', '1', '--parallel', '2', '--timeout', str(timeout)],
                              capture_output=True, text=True)

    def summary(self):
        return json.loads((self.output / 'batch_1/run_summary.json').read_text())

    def test_success_preserves_finite_approximation_errors(self):
        result = self.launch(['ok', 'range'])
        self.assertEqual(result.returncode, 0, result.stderr)
        summary = self.summary()
        self.assertTrue(summary['completed'])
        self.assertGreater(summary['records'][1]['outside_tolerance_1e_12'], 0)
        m = json.loads((self.output / 'batch_1/manifest.json').read_text())
        self.assertEqual(m['binary']['build_info']['source_revision'], 'a'*40)
        self.assertTrue((self.output / 'batch_1/.inputs/LinRacc').is_file())

    def test_backend_failure_propagates_through_shell_job(self):
        result = self.launch(['ok', 'exit_fail'], wrapper=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.summary()['success'], 1)
        self.assertEqual(self.summary()['failed'], 1)

    def test_zero_exit_without_output_is_failure(self):
        self.assertNotEqual(self.launch(['missing']).returncode, 0)
        self.assertFalse(self.summary()['completed'])

    def test_truncated_output_is_failure(self):
        self.assertNotEqual(self.launch(['truncated']).returncode, 0)
        self.assertFalse((self.output / 'batch_1/truncated.txt').exists())

    def test_nonfinite_output_is_failure(self):
        self.assertNotEqual(self.launch(['nan']).returncode, 0)
        self.assertEqual(self.summary()['failed'], 1)

    def test_revision_mismatch_is_failure(self):
        self.assertNotEqual(self.launch(['wrong_revision']).returncode, 0)
        self.assertIn('metadata mismatch', self.summary()['records'][0]['error'])

    def test_timeout_is_recorded_and_failed(self):
        self.assertNotEqual(self.launch(['timeout'], timeout=.1).returncode, 0)
        self.assertEqual(self.summary()['records'][0]['status'], 'TIMEOUT')

    def test_duplicate_ids_fail_before_processing(self):
        self.assertNotEqual(self.launch(['ok', 'ok']).returncode, 0)
        self.assertFalse((self.output / 'batch_1').exists())

    def test_filename_collisions_fail_before_processing(self):
        self.assertNotEqual(self.launch(['a/b', 'a_b']).returncode, 0)
        self.assertFalse((self.output / 'batch_1').exists())

    def test_empty_chunk_is_not_success(self):
        self.assertNotEqual(self.launch([]).returncode, 0)

    def test_dirty_build_is_rejected(self):
        self.binary.write_text(FAKE.replace("revision = 'a' * 40", "revision = 'a' * 40 + '-dirty'"))
        result = self.launch(['ok'])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('clean Git revision', result.stderr)

    def test_unknown_build_is_rejected(self):
        self.binary.write_text(FAKE.replace("revision = 'a' * 40", "revision = 'unknown'"))
        result = self.launch(['ok'])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('clean Git revision', result.stderr)

    def test_existing_results_are_not_overwritten(self):
        first = self.launch(['ok'])
        self.assertEqual(first.returncode, 0, first.stderr)
        before = (self.output / 'batch_1/manifest.json').read_bytes()
        self.assertNotEqual(self.launch(['ok']).returncode, 0)
        self.assertEqual(before, (self.output / 'batch_1/manifest.json').read_bytes())


if __name__ == '__main__':
    unittest.main()
