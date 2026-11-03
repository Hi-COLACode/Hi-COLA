"""Analytic and finite-difference checks of the numerical frontend."""
import unittest
import numpy as np
from HiCOLA.Frontend import numerical_solver as ns


class CosmologyTests(unittest.TestCase):
    def test_expansion_and_density_reference(self):
        # Pure component universes have exact power-law expansion histories.
        for dtype in (int, float, np.float64):
            for r, m, power in ((1, 0, 2), (0, 1, 1.5), (0, 0, 0)):
                for z in (dtype(0), dtype(1), np.array([0., 1., 3.])):
                    expected = (1 + z) ** power
                    for fn in (ns.comp_E_LCDM, ns.comp_E_LCDM_DE):
                        np.testing.assert_allclose(fn(z, r, m), expected)
                    np.testing.assert_allclose(ns.comp_H_LCDM(z, r, m, 70.), 70 * expected)
                    densities = [ns.comp_Omega_r_LCDM(z, r, m),
                                 ns.comp_Omega_m_LCDM(z, r, m),
                                 ns.comp_Omega_L_LCDM(z, r, m)]
                    for actual, reference in zip(densities, (r, m, 1-r-m)):
                        np.testing.assert_allclose(actual, reference)

    def test_mixed_universe_normalization(self):
        z = np.array([0., 0.2, 1., 10.])
        densities = [fn(z, 0.001, 0.3) for fn in
                     (ns.comp_Omega_r_LCDM, ns.comp_Omega_m_LCDM, ns.comp_Omega_L_LCDM)]
        np.testing.assert_allclose(sum(densities), 1.)
        self.assertEqual(ns.comp_E_LCDM(0., 0.001, 0.3), 1.)
        np.testing.assert_allclose(ns.comp_Omega_DE_LCDM(-np.log1p(z), .001, .3), densities[2])

    def test_derivatives_against_finite_differences(self):
        x = np.array([-2., -1., 0.])
        dx = 1e-5
        r, m = .001, .3
        E = ns.comp_E_LCDM(np.exp(-x)-1, r, m)
        EpE = ns.comp_E_prime_E_LCDM(x, r, m)
        logE = lambda t: np.log(ns.comp_E_LCDM(np.exp(-t)-1, r, m))
        np.testing.assert_allclose(EpE, (logE(x+dx)-logE(x-dx))/(2*dx), rtol=1e-8)
        for density, derivative in ((ns.comp_Omega_r_LCDM, ns.comp_Omega_r_prime),
                                    (ns.comp_Omega_m_LCDM, ns.comp_Omega_m_prime)):
            reference = (density(np.exp(-x-dx)-1,r,m)-density(np.exp(-x+dx)-1,r,m))/(2*dx)
            np.testing.assert_allclose(derivative(density(np.exp(-x)-1,r,m), E, E*EpE), reference, rtol=1e-8)
        reference = (ns.comp_Omega_DE_LCDM(x+dx,r,m)-ns.comp_Omega_DE_LCDM(x-dx,r,m))/(2*dx)
        np.testing.assert_allclose(ns.comp_Omega_DE_prime_LCDM(EpE, ns.comp_Omega_DE_LCDM(x,r,m)), reference, rtol=1e-8)
        np.testing.assert_allclose(ns.comp_Omega_l_prime(1-r-m,E,E*EpE), reference, rtol=1e-8)

    def test_alpha_and_mass_integration(self):
        x = np.array([-1., 0.])
        de = ns.comp_Omega_DE_LCDM(x, 0., 0.)
        np.testing.assert_allclose(ns.comp_alpha_M_propto_Omega_DE_LCDM(.2, de), .2)
        np.testing.assert_allclose(ns.alpha_M_int_propto_Omega_DE_LCDM(x,0.,0.,.2), .2)
        np.testing.assert_allclose(ns.comp_alpha_M_prime_propto_Omega_DE_LCDM(.2,np.zeros(2)), 0.)
        # Constant alpha in a pure Lambda universe integrates to c_M * interval.
        self.assertAlmostEqual(ns.comp_Meffsq_x2_x1_propto_Omega_DE_LCDM(-2.,0.,0.,0.,.2), np.exp(.4))
        self.assertAlmostEqual(ns.comp_Meffsq_x2_x1_propto_Omega_DE_LCDM(0.,-2.,0.,0.,.2), np.exp(-.4))

    def test_invalid_redshift_type(self):
        for fn in (ns.comp_E_LCDM, ns.comp_E_LCDM_DE, ns.comp_Omega_r_LCDM,
                   ns.comp_Omega_m_LCDM, ns.comp_Omega_L_LCDM):
            with self.subTest(function=fn.__name__), self.assertRaises(TypeError):
                fn(object(), .001, .3)

    def test_screening_accepts_lists_and_arrays(self):
        for container in (list, np.array):
            np.testing.assert_allclose(ns.chi_over_delta(container([1.,.5]), container([1.,2.]),
                container([2.,2.]), container([3.,3.]), .3), [1.8,3.6])
        with self.assertRaises(ValueError):
            ns.chi_over_delta([1,2,3], [1,2], [1], [1], .3)

    def test_friedmann_closure_reference(self):
        # No scalar contribution: normalized matter plus radiation closes exactly.
        self.assertEqual(ns.fried1(0., 2., 3., .2, .8, 1., 0., 1., 1.), 0.)


class SolverTests(unittest.TestCase):
    def test_closure_argument_substitution(self):
        for index in range(0,5):
            with self.subTest(index=index):
                result = ns.fried_RHS_wrapper(9., ['odeint_parameters',index],
                    lambda *args: args, 1.,2.,3.,4.,5.,6.,[7.])
                expected = [1.,2.,3.,4.,5.,6.,7.]
                expected[index] = 9.
                self.assertEqual(result, tuple(expected))
        self.assertEqual(ns.fried_RHS_wrapper(9., ['parameters',0],lambda *a:a,
                         1.,2.,3.,4.,5.,6.,[7.]), (1.,2.,3.,4.,5.,6.,9.))

    def test_invalid_closure_declaration(self):
        for declaration in (['unknown',0], ['odeint_parameters',99]):
            with self.assertRaisesRegex(Exception, 'Invalid'):
                ns.fried_RHS_wrapper(1.,declaration,lambda *a:0.,1.,2.,3.,4.,5.,6.,[])

    def test_root_solver_against_exact_linear_root(self):
        for index in range(0,5):
            def residual(*args):
                return args[index] - .75
            root, diagnostics = ns.comp_param_close(residual,['odeint_parameters',index],
                                   1.,0.,.1,.01,.3,.69,[])
            self.assertAlmostEqual(root, .75)
            self.assertEqual(diagnostics[2], 1)

    def test_primes_normal_safe_and_gr_paths(self):
        for A in (1., 1e-5, -1e-5):
            calls = []
            def normal(*args):
                calls.append('normal')
                return -1.5
            def safe(*args):
                calls.append(('safe',args[-1]))
                return -1.5
            result = ns.comp_primes(0., [0.,0.,1.,0.,1.,0.], 1.,0.,1.,0.,
                normal,safe,lambda *a:0.,lambda *a:0.,lambda *a:A,['parameters',0],[])
            np.testing.assert_allclose(result,[0.,0.,-1.5,0.,0.,0.])
            self.assertEqual(calls, ['normal'] if abs(A)>=1e-3 else [('safe',np.sign(A))])

    def test_run_solver_matter_universe(self):
        # Exercise SciPy integration through the real wrapper; EdS has E=a^-3/2.
        data = dict(cosmological_parameters=[0.,1.,0.], initial_conditions=[1.,0.,0.],
            simulation_parameters=[30,3.,False,0.,True], Horndeski_parameters=[],
            GG4_over_GN=lambda *a:1., EpE_debug={}, omega_phi_debug={},
            Horndeski_K=lambda *a:0., forwards_flag=False,
            symbol_list=[], odeint_parameter_symbols=[], closure_declaration=['odeint_parameters',1],
            fried_RHS_lambda=lambda E,phi,p,r,m,l:phi,
            E_prime_E_lambda=lambda *a:-1.5, E_prime_E_safelambda=lambda *a:-1.5,
            phi_primeprime_lambda=lambda *a:0., phi_primeprime_safelambda=lambda *a:0.,
            omega_phi_lambda=lambda *a:0., A_lambda=lambda *a:1.,
            calB_lambda=lambda *a:0., calC_lambda=lambda *a:0., coupling_factor=lambda *a:0.)
        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=UserWarning, module=r'^scipy(\.|$)')
            result = ns.run_solver(data)
        np.testing.assert_allclose(result['Hubble'], result['a']**-1.5, rtol=2e-6)
        np.testing.assert_allclose(result['omega_m'],1.)
        for name in ('scalar_prime','scalar_primeprime','chi_over_delta','coupling_factor'):
            np.testing.assert_allclose(result[name],0.)
        self.assertEqual(result['a'][0],1.)
        self.assertAlmostEqual(result['a'][-1],.25)
