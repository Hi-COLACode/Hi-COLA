"""Configuration conversion tests using temporary copies of shipped examples."""
import tempfile
import unittest
from pathlib import Path
from configobj import ConfigObj
import numpy as np
from HiCOLA.Frontend import read_parameters as rp

ROOT = Path(__file__).resolve().parents[1] / 'HiCOLA'


def scan_settings():
    return ConfigObj(dict(model_name='test', cosmo_name='test',
        read_scan_parameters_from_file='False', suppression_flag='False',
        max_redshift='3', Npoints='30', GR_flag='True', threshold='0',
        early_DE_threshold='0.1', yellow_switch='True', blue_switch='True',
        red_switch='True', tolerance='0.01', Omega_m_crit='0.3',
        set_all_to_one='True', K='X', G3='0', G4='0.5',
        closure_guess_value='0', use_constraint_eq_on_odeint_variables='True',
        which_odeint_par='1'))


class ParameterTests(unittest.TestCase):
    def test_scan_settings_types(self):
        with tempfile.TemporaryDirectory() as directory:
            config = scan_settings()
            config.filename = str(Path(directory)/'settings.ini')
            config.write()
            result = rp.read_in_scan_settings(config.filename)
        self.assertIsInstance(result['simulation_parameters'][0],int)
        self.assertIsInstance(result['simulation_parameters'][1],float)
        self.assertIsInstance(result['simulation_parameters'][2],bool)
        self.assertEqual(result['mass_ratio_list'],[1.]*7)

    def test_shipped_parameters(self):
        result = rp.read_in_parameters(str(ROOT/'Parameter_files/Frontend_parameter_files/examples/ESS/ESS-A.ini'),
                                       str(ROOT/'Parameter_files/Frontend_parameter_files/numerical_parameters.ini'))
        r,m,l = result['cosmological_parameters']
        self.assertAlmostEqual(r,4.28e-5/.7307**2)
        self.assertAlmostEqual(m,(.02196+.1274)/.7307**2)
        self.assertTrue(np.isfinite(l))
        self.assertIsInstance(result['simulation_parameters'][0],int)

    def test_missing_keys_and_invalid_values(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'settings.ini'
            path.write_text('model_name = test\n')
            with self.assertRaises(KeyError):
                rp.read_in_scan_settings(str(path))
            for key,value in (('Npoints','invalid'),('suppression_flag','invalid')):
                config = scan_settings()
                config.filename = str(path)
                config[key] = value
                config.write()
                with self.subTest(key=key), self.assertRaises(ValueError):
                    rp.read_in_scan_settings(str(path))
