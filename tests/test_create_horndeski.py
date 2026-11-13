"""Checks of the numerical functions returned by expression_builder.create_Horndeski."""
import unittest

import sympy as sym

from HiCOLA.Frontend import expression_builder as eb

k1, g31 = sym.symbols('k_1 g_31')


class CreateHorndeskiLambdaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Cubic Galileon: K = k1 X, G3 = g31 X, G4 = 1/2 with unit mass ratios
        cls.fns = eb.create_Horndeski(k1*eb.X, g31*eb.X, sym.Rational(1, 2),
                                      [k1, g31], [1]*7)

    def test_horndeski_function_lambdas(self):
        # Arguments are (E, phi, phiprime, k1, g31); X = E^2 phiprime^2 / 2 = 0.5
        args = (2.0, 0.0, 0.5, 3.0, 5.0)
        self.assertAlmostEqual(self.fns['Horndeski_K'](*args), 1.5)
        self.assertAlmostEqual(self.fns['Horndeski_G3'](*args), 2.5)
        self.assertAlmostEqual(self.fns['Horndeski_G4'](*args), 0.5)

    def test_beta0_lambda_is_numeric(self):
        value = self.fns['beta0_lambda'](1.2, -0.3, 0.1, 0.4, 0.05, -1.0, 0.5)
        self.assertIsInstance(float(value), float)


if __name__ == '__main__':
    unittest.main()
