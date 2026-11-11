"""Check the regression-only figures and tolerance units."""
from pathlib import Path
import tempfile
import unittest

import numpy as np
from plots import overview_plot, tolerance_ratio, velocity_plot, worst_ratio


def is_png(path):
    return Path(path).read_bytes().startswith(b'\x89PNG')


class PlotTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)

    def test_ratio_is_in_tolerance_units(self):
        reference = np.array([1., 100.])
        # 'power' allows 1e-6 + 1e-4*|baseline|, so these are exactly one allowance away.
        actual = reference+1e-6+1e-4*reference
        np.testing.assert_allclose(tolerance_ratio(actual, reference, 'power'), [1., 1.])
        self.assertAlmostEqual(worst_ratio(actual, reference, 'power'), 1.)
        with self.assertRaises(ValueError):
            worst_ratio(actual[:1], reference, 'power')

    def test_velocity_figure_within_and_beyond_tolerance(self):
        stats = np.array([0., 0., 0., 400.])
        for name, actual in (('same', stats), ('drift', stats*[1, 1, 1, 1.01])):
            with self.subTest(name=name):
                self.assertTrue(is_png(velocity_plot(actual, stats, 'gr-np1', self.root/f'{name}.png')))

    def test_overview_figure_with_exact_and_failing_checks(self):
        path = overview_plot({'gr-np1 / power spectra': 0., 'gr-np1 / growth tables': 3.},
                             self.root/'overview.png')
        self.assertTrue(is_png(path))


if __name__ == '__main__':
    unittest.main()
