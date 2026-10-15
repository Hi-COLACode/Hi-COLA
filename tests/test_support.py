import os
import tempfile
import unittest
from pathlib import Path
import numpy as np
from configobj import ConfigObj
from HiCOLA.Utilities.Other import support as s, suppressor


class SupportTests(unittest.TestCase):
    def test_scan_arrays(self):
        for start,end,count,expected in ((0,1,3,[0,.5,1]),(2,2,9,[2]),(1,0,3,[1,.5,0]),(0,1,0,[])):
            np.testing.assert_allclose(s.make_scan_array(start,end,count),expected)
        config = ConfigObj({'q_min':'0','q_max':'1','q_number':'3'})
        np.testing.assert_allclose(s.generate_scan_array(config,'q'),[0,.5,1])
        for count,error in ((-1,ValueError),(1.5,TypeError)):
            with self.assertRaises(error):
                s.make_scan_array(0.,1.,count)
        with self.assertRaises(KeyError):
            s.generate_scan_array(config,'missing')
        config['q_number']='invalid'
        with self.assertRaises(ValueError):
            s.generate_scan_array(config,'q')

    def test_nearest_index(self):
        for arr in ([1.,3.,5.],np.array([1.,3.,5.])):
            self.assertEqual(s.nearest_index(arr,3.2),1)
            self.assertEqual(s.nearest_index(arr,2.),0)
            self.assertEqual(s.nearest_index(arr,100.),2)
        with self.assertRaises(ValueError):
            s.nearest_index([],1.)
        with self.assertRaises(TypeError):
            s.nearest_index([1.,2.],object())

    def test_output_tables(self):
        writers = ((s.write_data,6),(s.write_data_coupl,5),(s.write_data_screencoupl,3),(s.write_all_data,11))
        with tempfile.TemporaryDirectory() as directory:
            for fn,count in writers:
                columns = [[i+.123456, i+1.234567] for i in range(count)]
                path = Path(directory)/'table.txt'
                fn(*columns,path)
                np.testing.assert_allclose(np.loadtxt(path),np.array(columns).T[::-1],rtol=5e-5)
            flex_columns = [[.123456, 1.234567], [2.345678, 3.456789]]
            flex_path = Path(directory)/'flex_table.txt'
            s.write_data_flex(flex_columns,flex_path)
            np.testing.assert_allclose(np.loadtxt(flex_path),np.array(flex_columns).T[::-1],rtol=5e-5)

    def test_output_missing_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(FileNotFoundError):
                s.write_data_flex([[1.,2.]],Path(directory)/'missing'/'table.txt')

    def test_ess_parameters_constraints(self):
        for fraction in (0.,.5,1.):
            U,r,m,l,(k1,k2,g31,g32) = s.ESS_dS_parameters(2.,fraction,-1.,3.,h=1.,Omega_r0h2=.1,Omega_b0h2=.1,Omega_c0h2=.2)
            self.assertEqual(U,.5)
            self.assertAlmostEqual(r+m+l,1.-fraction*.6)
            alpha = 1-l/4
            self.assertAlmostEqual(k1,-alpha)
            self.assertAlmostEqual(g31,3*alpha)
            self.assertAlmostEqual(k2+2*k1,-12*alpha)
            self.assertAlmostEqual(l/4+k1/6+k2/4+g31+2*g32,1.)

    def test_renamer(self):
        self.assertEqual(s.renamer('sample.txt'),'sample_1.txt')
        self.assertEqual(s.renamer('sample_1.txt'),'sample_2.txt')

    def test_fileno(self):
        self.assertEqual(suppressor.fileno(1),1)
        with tempfile.TemporaryFile() as file:
            self.assertEqual(suppressor.fileno(file),file.fileno())
        for value in ('abc',None,1.5):
            with self.assertRaisesRegex(ValueError,'file descriptor'):
                suppressor.fileno(value)

    def test_descriptor_redirect_restored_after_exception(self):
        with tempfile.TemporaryFile(mode='w+b') as source, tempfile.TemporaryFile(mode='w+b') as target:
            with self.assertRaisesRegex(RuntimeError,'test'):
                with suppressor.stdout_redirected(target,stdout=source):
                    os.write(source.fileno(),b'redirected')
                    raise RuntimeError('test')
            os.write(source.fileno(),b'restored')
            source.seek(0)
            target.seek(0)
            self.assertEqual(source.read(),b'restored')
            self.assertEqual(target.read(),b'redirected')
