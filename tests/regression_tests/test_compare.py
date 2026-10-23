"""Check that the regression gate rejects broken artifacts and numerical drift."""
import io
from pathlib import Path
import struct
import tempfile
import unittest
import numpy as np
from compare import close, force, lightcone, numeric, particles, record


class ComparisonTests(unittest.TestCase):
    def test_tolerance_boundary(self):
        ref = np.array([1.,0.])
        close(ref+np.array([1e-6,1e-9]),ref,'tables')
        with self.assertRaises(AssertionError):
            close(ref+.1,ref,'tables')

    def test_rejects_nonfinite_empty_and_shape_mismatch(self):
        for actual in (np.array([np.nan]),np.array([np.inf]),np.array([]),np.ones(2)):
            with self.subTest(actual=actual), self.assertRaises(ValueError):
                close(actual,np.ones(1),'tables')
        with self.assertRaises(ValueError):
            close(np.ones(1),np.array([np.nan]),'tables')

    def test_numeric_rejects_nan_and_wrong_columns(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'data.txt'
            path.write_text('1 nan\n')
            with self.assertRaises(ValueError):
                numeric(path)
            path.write_text('1 2\n')
            with self.assertRaises(ValueError):
                numeric(path,3)

    def test_force_requires_complete_unique_grid(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'force.txt'
            np.savetxt(path,np.zeros((16**3,6)))
            with self.assertRaisesRegex(ValueError,'Duplicate'):
                force([path])
            np.savetxt(path,np.zeros((1,6)))
            with self.assertRaisesRegex(ValueError,'Incomplete'):
                force([path])

    def test_lightcone_requires_nonempty_cone(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'cone.txt'
            path.write_text('# empty rank\n')
            with self.assertRaisesRegex(ValueError,'No lightcone'):
                lightcone([path])
            path.write_text('3 2 1\n1 2 3\n')
            np.testing.assert_array_equal(lightcone([path]),[[1,2,3],[3,2,1]])

    def test_record_markers_and_truncation(self):
        for raw in (b'',struct.pack('<I',5)+b'bad',struct.pack('<I',1)+b'x'+struct.pack('<I',2)):
            with self.assertRaises(ValueError):
                record(io.BytesIO(raw))

    def test_snapshot_epoch_and_particle_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'gadget'
            n = 16**3
            header = bytearray(256)
            struct.pack_into('<I',header,4,n)
            struct.pack_into('<2d',header,72,1.,0.)
            def write(time):
                struct.pack_into('<d',header,72,time)
                with path.open('wb') as stream:
                    for data in (header,np.zeros((n,3),dtype='<f4').tobytes(),
                                 np.ones((n,3),dtype='<f4').tobytes(),np.arange(n,dtype='<u8')[::-1].tobytes()):
                        size = struct.pack('<I',len(data))
                        stream.write(size+data+size)
            write(1.)
            ids,pos,vel = particles([path])
            np.testing.assert_array_equal(ids,np.arange(n))
            self.assertEqual(pos.shape,(n,3))
            np.testing.assert_array_equal(vel,np.ones((n,3)))
            write(.5)
            with self.assertRaisesRegex(ValueError,'z=0'):
                particles([path])


if __name__=='__main__':
    unittest.main()
