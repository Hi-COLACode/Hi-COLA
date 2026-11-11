"""Run the hicola_plot_* tools on small synthetic Hi-COLA output folders."""
import contextlib
import io
from pathlib import Path
import tempfile
import unittest

import numpy as np

from HiCOLA.Plotting import background, frontend, growth, pofk
from HiCOLA.Plotting.common import read_table


def write(path, header, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, np.column_stack(columns), header=header)


def make_output(folder, name, scale=1.):
    """Write the backend tables the tools read, for redshifts 0 and 1 and two k values."""
    k = np.logspace(-1, 0, 8)
    for z in (0., 1.):
        write(folder/f'pofk_{name}_cb_z{z:.3f}.txt',
              'k  (h/Mpc)    Pcb(k)  (Mpc/h)^3    Pcb_linear(k) (Mpc/h)^3    Pcb_NL(k) (Mpc/h)^3   '
              'ShotnoiseSubtracted = true', (k, scale*1e3/k/(1+z), 1e3/k, 1e3/k))
    a = np.geomspace(1e-4, 10, 50)
    H = np.sqrt(.3/a**3+.7)
    write(folder/f'cosmology_{name}.txt', 'a H/H0 dlogH/dloga OmegaM OmegaR OmegaNu OmegaMNu OmegaNu_exact OmegaLambda',
          (a, scale*H, -1.5+0*a, .3/a**3/H**2, 0*a, 0*a, 0*a, 0*a, .7/H**2))
    for wavenumber in ('0.1', '1.0'):
        write(folder/f'gravitymodel_{name}_k{wavenumber}.txt', 'a GeffG(a,k) GG4_GN(a) D1(a,k)',
              (a, 1+0*a, 1+0*a, scale*a))


def make_frontend(folder, name, force_columns=4):
    a = np.geomspace(1e-4, 1, 40)
    folder.mkdir(parents=True)
    np.savetxt(folder/f'{name}_test_expansion.txt', np.c_[a, np.sqrt(.3/a**3+.7), -1.5+0*a])
    np.savetxt(folder/f'{name}_test_force.txt', np.c_[a, a, .1*a, 1+0*a][:, :force_columns])


def run(tool, *argv):
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        return tool.main([str(value) for value in argv])


class PlottingToolTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        make_output(self.root/'gr'/'output', 'gr')
        make_output(self.root/'cubic'/'output', 'cubic', scale=1.02)

    def assert_png(self, path):
        self.assertTrue(path.read_bytes().startswith(b'\x89PNG'), path)

    def test_single_run_figures_are_written_beside_the_output_folder(self):
        for tool, name in ((pofk, 'pofk'), (background, 'background'), (growth, 'growth')):
            with self.subTest(tool=name):
                self.assertEqual(run(tool, self.root/'gr'/'output'), 0)
                self.assert_png(self.root/'gr'/'plots'/f'{name}.png')

    def test_reference_runs_add_ratio_panels(self):
        for tool, name in ((pofk, 'pofk'), (background, 'background'), (growth, 'growth')):
            with self.subTest(tool=name):
                run(tool, self.root/'cubic'/'output', '--reference', self.root/'gr'/'output',
                    '--bands', '0.01', '--out', self.root/'boost')
                self.assert_png(self.root/'boost'/f'{name}.png')

    def test_frontend_tables_with_three_or_four_force_columns(self):
        make_frontend(self.root/'frontend'/'cubic', 'cubic')
        make_frontend(self.root/'frontend'/'kmouflage', 'kmouflage', force_columns=3)
        run(frontend, self.root/'frontend'/'kmouflage')
        run(frontend, self.root/'frontend'/'cubic', '--reference', self.root/'frontend'/'kmouflage',
            '--format', 'pdf')
        self.assert_png(self.root/'frontend'/'kmouflage'/'plots'/'frontend.png')
        self.assertTrue((self.root/'frontend'/'cubic'/'plots'/'frontend.pdf').is_file())

    def test_at_most_two_bands(self):
        for bands in (['0.01', '0.02', '0.03'], ['-0.01']):
            with self.subTest(bands=bands), self.assertRaises(SystemExit):
                run(pofk, self.root/'gr'/'output', '--bands', *bands)
        run(pofk, self.root/'cubic'/'output', '--reference', self.root/'gr'/'output', '--bands')
        self.assert_png(self.root/'cubic'/'plots'/'pofk.png')

    def test_missing_inputs_are_reported(self):
        (self.root/'empty').mkdir()
        for tool in (pofk, background, growth):
            with self.subTest(tool=tool.__name__), self.assertRaises(SystemExit):
                run(tool, self.root/'empty')
        with self.assertRaises(SystemExit):
            run(pofk, self.root/'cubic'/'output', '--reference', self.root/'empty')

    def test_redshift_selection(self):
        self.assertEqual(list(pofk.spectra(self.root/'gr'/'output', [0.])), ['0.000'])
        self.assertEqual(list(pofk.spectra(self.root/'gr'/'output')), ['1.000', '0.000'])

    def test_header_names_keep_units(self):
        table = read_table(self.root/'gr'/'output'/'pofk_gr_cb_z0.000.txt')
        self.assertEqual(list(table), ['k (h/Mpc)', 'Pcb(k) (Mpc/h)^3',
                                       'Pcb_linear(k) (Mpc/h)^3', 'Pcb_NL(k) (Mpc/h)^3'])


if __name__ == '__main__':
    unittest.main()
