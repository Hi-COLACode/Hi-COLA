"""Check that summary outputs are measured against their tolerances from a synthetic run folder."""
from pathlib import Path
import struct
import tempfile
import unittest

import numpy as np
from collect import summarise


def write_gadget(path, velocities):
    def record(payload):
        size = struct.pack('<I', len(payload))
        return size+payload+size
    count = len(velocities)
    header = bytearray(256)
    struct.pack_into('<6I', header, 0, count, 0, 0, 0, 0, 0)
    struct.pack_into('<2d', header, 72, 1.0, 0.0)
    with open(path, 'wb') as stream:
        stream.write(record(bytes(header)))
        stream.write(record(np.zeros((count, 3), '<f4').tobytes()))
        stream.write(record(velocities.astype('<f4').tobytes()))
        stream.write(record(np.arange(count, dtype='<u4').tobytes()))


class CollectTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.work = Path(directory.name)

    def make_run(self, power_scale=1.01, candidate_rows=8):
        k = np.logspace(-2, 0, 8)
        for label in ('baseline', 'candidate'):
            output = self.work/label/'gr-np1'/'output'
            (output/'snapshot_gr_z0.000').mkdir(parents=True)
            scale = power_scale if label == 'candidate' else 1.
            rows = candidate_rows if label == 'candidate' else 8
            np.savetxt(output/'pofk_gr_cb_z0.000.txt', np.c_[k, scale/k][:rows])
            np.savetxt(output/'snapshot_gr_z0.000'/'pofk_z0.000.txt', np.c_[k, 1/k])
            np.savetxt(output/'gravitymodel_gr_k0.1.txt', np.c_[k, k])
            frontend = self.work/label/'frontend'/'kmouflage'
            frontend.mkdir(parents=True)
            np.savetxt(frontend/'kmouflage_test_expansion.txt', np.c_[k, k, k])
            write_gadget(output/'snapshot_gr_z0.000'/'gadget_z0.000.0',
                         np.tile([3.0, 0.0, 4.0], (16**3, 1)))

    def test_every_output_is_measured(self):
        self.make_run()
        result = summarise(self.work, ['gr'], [1], self.work/'plots')
        self.assertEqual(result['errors'], {})
        written = sorted(str(p.relative_to(self.work/'plots')) for p in (self.work/'plots').rglob('*.png'))
        self.assertEqual(written, ['gr-np1/velocity_statistics.png', 'overview.png'])
        self.assertEqual(result['checks']['frontend / kmouflage'], 0)
        # The 1% power difference is 100 tolerance units.
        self.assertGreater(result['checks']['gr-np1 / power spectra'], 1)
        self.assertEqual(result['checks']['gr-np1 / growth tables'], 0)
        self.assertEqual(result['checks']['gr-np1 / velocity statistics'], 0)

    def test_problems_are_recorded_and_other_checks_continue(self):
        self.make_run(candidate_rows=5)
        result = summarise(self.work, ['gr', 'kmouflage'], [1], self.work/'plots')
        self.assertIn('gr-np1/pofk_gr_cb_z0.000.txt', result['errors'])
        self.assertIn('kmouflage-np1/velocity_statistics', result['errors'])
        self.assertIn('gr-np1 / growth tables', result['checks'])
        self.assertTrue((self.work/'plots'/'overview.png').is_file())


if __name__ == '__main__':
    unittest.main()
