"""Hand-derived polynomial references for symbolic derivative helpers."""
import unittest
import sympy as sym
from HiCOLA.Frontend import expression_builder as eb


class ExpressionTests(unittest.TestCase):
    def test_polynomial_derivatives(self):
        X, phi = eb.X, eb.phi
        polynomial = X**3 * phi**2 + 2*phi
        expected = (3*X**2*phi**2,6*X*phi**2,6*X**2*phi,6*X**2*phi,2*X**3,2*X**3*phi+2)
        for fn in (eb.K_func,eb.G3_func,eb.G4_func):
            for selector, reference in zip(('X','XX','X,phi'), expected[:3]):
                with self.subTest(function=fn.__name__, selector=selector):
                    self.assertEqual(sym.expand(fn(polynomial,selector)-reference),0)
            references = (expected[0],expected[1],expected[2],expected[5]) if fn==eb.K_func else expected
            self.assertEqual(fn(polynomial), references)
        self.assertEqual(eb.K_func(polynomial,'phi'),expected[-1])
        for fn in (eb.G3_func,eb.G4_func):
            self.assertEqual(fn(polynomial,'phi,X'),expected[3])
            self.assertEqual(fn(polynomial,'phi,phi'),expected[4])

    def test_constants_and_symbols(self):
        for fn in (eb.K_func,eb.G3_func,eb.G4_func):
            for value in (1,1.,sym.Integer(1)):
                self.assertTrue(all(v==0 for v in fn(value)))
        self.assertEqual(eb.sy('q'),sym.Symbol('q'))
        namespace = {'sym':sym}
        exec(eb.declare_symbols(),namespace)
        self.assertEqual(namespace['X'],eb.X)

    def test_invalid_expression_types(self):
        for fn in (eb.K_func,eb.G3_func,eb.G4_func,eb.sy):
            with self.subTest(function=fn.__name__), self.assertRaises((TypeError,sym.SympifyError)):
                fn(object())

    @unittest.expectedFailure
    def test_g3_phi_derivative(self):
        self.assertEqual(eb.G3_func(eb.X*eb.phi**2,'phi'),2*eb.X*eb.phi)
