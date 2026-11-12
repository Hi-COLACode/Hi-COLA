"""Exercise summary and particle checks against independent output fixtures."""
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np
from run import compare_runs


class FlavourTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.baseline = self.root/'baseline'
        self.candidate = self.root/'candidate'
        self.args = SimpleNamespace(cases=['gr'], ranks=[1])
        self.positions = np.arange(16**3*3, dtype=np.float32).reshape(-1, 3)/4096
        self.velocities = np.tile([[3., 4., 0.], [-3., -4., 0.]], (16**3//2, 1))
        for root in (self.baseline, self.candidate):
            self.make_case(root, 'gr')

    def make_case(self, root, case):
        folder = root/f'{case}-np1'
        output = folder/'output'
        output.mkdir(parents=True)
        np.savetxt(output/'gravitymodel_gr.txt', [[0.1, 1.], [1., 1.]])
        np.savetxt(output/'pofk_gr.txt', [[0.1, 10.], [0.2, 20.]])
        self.snapshot(output/'gadget_z0.000.0', self.positions, self.velocities)
        coords = np.indices((16, 16, 16)).reshape(3, -1).T/16
        np.savetxt(folder/'force-False.0.txt', np.column_stack((coords, np.zeros_like(coords))))

    def snapshot(self, path, positions, velocities):
        count = len(positions)
        header = bytearray(256)
        struct.pack_into('<I', header, 4, count)
        struct.pack_into('<2d', header, 72, 1., 0.)
        with path.open('wb') as stream:
            for data in (header, np.asarray(positions, dtype='<f4').tobytes(),
                         np.asarray(velocities, dtype='<f4').tobytes(),
                         np.arange(count, dtype='<u8').tobytes()):
                size = struct.pack('<I', len(data))
                stream.write(size+data+size)

    def compare(self, flavour):
        return compare_runs(self.baseline, self.candidate, self.args, flavour)

    def test_unchanged_outputs_pass_both_flavours(self):
        for flavour in ('summary', 'particles'):
            with self.subTest(flavour=flavour):
                self.assertFalse(self.compare(flavour)['failures'])

    def test_particle_translation_does_not_fail_summary(self):
        self.snapshot(self.candidate/'gr-np1/output/gadget_z0.000.0',
                      self.positions+1., self.velocities)
        self.assertFalse(self.compare('summary')['failures'])
        self.assertIn('gr-np1/positions', self.compare('particles')['failures'])

    def test_reassigning_velocities_preserves_statistics_but_changes_particles(self):
        self.snapshot(self.candidate/'gr-np1/output/gadget_z0.000.0',
                      self.positions, np.roll(self.velocities, 1, axis=0))
        self.assertFalse(self.compare('summary')['failures'])
        self.assertIn('gr-np1/velocities', self.compare('particles')['failures'])

    def test_velocity_dispersion_drift_is_a_summary_failure(self):
        self.snapshot(self.candidate/'gr-np1/output/gadget_z0.000.0',
                      self.positions, self.velocities*1.1)
        self.assertIn('gr-np1/velocity_statistics', self.compare('summary')['failures'])

    def test_power_drift_is_a_summary_failure(self):
        np.savetxt(self.candidate/'gr-np1/output/pofk_gr.txt', [[0.1, 11.], [0.2, 22.]])
        self.assertIn('gr-np1/pofk_gr.txt', self.compare('summary')['failures'])
        self.assertFalse(self.compare('particles')['failures'])

    def test_missing_power_is_a_summary_failure(self):
        (self.candidate/'gr-np1/output/pofk_gr.txt').unlink()
        self.assertIn('gr-np1/valid_outputs', self.compare('summary')['failures'])

    def test_nonfinite_snapshots_are_rejected_by_the_summary_gate(self):
        positions = self.positions.copy()
        positions[0, 0] = np.nan
        self.snapshot(self.candidate/'gr-np1/output/gadget_z0.000.0', positions, self.velocities)
        self.assertIn('gr-np1/valid_outputs', self.compare('summary')['failures'])

    def test_force_drift_is_particle_level(self):
        path = self.candidate/'gr-np1/force-False.0.txt'
        data = np.loadtxt(path)
        data[:, 3:] = 0.1
        np.savetxt(path, data)
        self.assertFalse(self.compare('summary')['failures'])
        self.assertIn('gr-np1/force-False', self.compare('particles')['failures'])

    def test_compare_only_reuses_outputs_and_reports_advisory_changes(self):
        self.snapshot(self.candidate/'gr-np1/output/gadget_z0.000.0',
                      self.positions+1., self.velocities)
        manifest = dict(simulations_complete=True, status='simulated', cases=['gr'], ranks=[1],
                        baseline=dict(commit='baseline'), candidate=dict(commit='candidate'))
        (self.root/'report.json').write_text(json.dumps(manifest))
        env = dict(os.environ, PATH='', GITHUB_ACTIONS='true', GITHUB_STEP_SUMMARY=str(self.root/'summary.md'))
        for flavour, code in [('summary', 0), ('particles', 1)]:
            with self.subTest(flavour=flavour):
                result = subprocess.run([sys.executable, str(Path(__file__).with_name('run.py')),
                    '--compare-only', '--work', str(self.root), '--flavour', flavour],
                    env=env, capture_output=True, text=True)
                self.assertEqual(result.returncode, code, result.stdout+result.stderr)
                report = json.loads((self.root/f'{flavour}-report.json').read_text())
                self.assertEqual(report['status'], 'passed' if code == 0 else 'failed')
                if code:
                    self.assertIn('::warning title=Particles regression differences::', result.stdout)
        self.assertEqual(json.loads((self.root/'report.json').read_text()), manifest)


if __name__ == '__main__':
    unittest.main()
