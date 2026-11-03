import unittest
import numpy as np
from compare import velocity_statistics


class TestVelocityStatistics(unittest.TestCase):
    def test_constant_vector(self):
        v = np.array([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]])
        expected = np.array([1.0, 2.0, 3.0, np.sqrt(1.0 + 4.0 + 9.0)])
        result = velocity_statistics(v)
        np.testing.assert_allclose(result, expected)

    def test_balanced_opposing_vectors(self):
        v = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
        expected = np.array([0.0, 0.0, 0.0, 1.0])
        result = velocity_statistics(v)
        np.testing.assert_allclose(result, expected)

    def test_reorder_invariance(self):
        v1 = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        v2 = np.array([[4.0, 5.0, 6.0], [1.0, 2.0, 3.0]])
        r1 = velocity_statistics(v1)
        r2 = velocity_statistics(v2)
        np.testing.assert_allclose(r1, r2)

    def test_empty_array_raises(self):
        v = np.empty((0, 3))
        with self.assertRaises(ValueError):
            velocity_statistics(v)

    def test_invalid_shapes_and_values_raise(self):
        wrong_shape_1d = np.array([1.0, 2.0, 3.0])
        with self.assertRaises(ValueError):
            velocity_statistics(wrong_shape_1d)

        wrong_shape_2d = np.array([[1.0, 2.0], [3.0, 4.0]])
        with self.assertRaises(ValueError):
            velocity_statistics(wrong_shape_2d)

        nan_array = np.array([[np.nan, 1.0, 2.0]])
        with self.assertRaises(ValueError):
            velocity_statistics(nan_array)

        inf_array = np.array([[1.0, np.inf, 2.0]])
        with self.assertRaises(ValueError):
            velocity_statistics(inf_array)


if __name__ == '__main__':
    unittest.main()
