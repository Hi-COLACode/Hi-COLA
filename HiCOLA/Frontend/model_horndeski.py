import numpy as np
import time
import sympy as sym
from scipy.optimize import newton

from .model_standard import StandardModel
from . import redshift, utils


class HorndeskiModel(StandardModel):

    """
    A class for constructing a user defined Horndeski gravity model and computing (numerically) the
    background expansion and growth functions to construct non-linear cosmological simulations using 
    the Hi-COLA Backend.
    """
    
    # Initialising the main class

    def __init__(self):
        """
        Initialises the Horndeski model class.
        """
        super().__init__()

        # Construct dictionary containing sympy symbolic variables.
        self.sym = {}
        
        # symbol for the scale factor
        self.sym['a'] = sym.symbols('a')

        # Main functions of time
        self.sym['E'] = sym.symbols('E')
        self.sym['phi'] = sym.symbols('phi')
        self.sym['phi_prime'] = sym.symbols("phi^{'}")
        self.sym['X'] = sym.symbols('X')

        # Renormalise constant for the normalised expansion rate E
        self.sym['f_H'] = sym.symbols('f_H')

        # Mass scalings
        self.sym['M_sp'] = sym.symbols('M_{sp}')
        self.sym['M_Kp2'] = sym.symbols('M_{Kp2}')
        self.sym['M_G3p'] = sym.symbols('M_{G3p}')
        self.sym['M_G4p2'] = sym.symbols('M_{G4p2}')

        # Photons
        self.sym['rho_g'] = sym.symbols('rho_g')
        self.sym['rhohat_g'] = sym.symbols('\hat{rho}_g')
        self.sym['Omega_g'] = sym.symbols('Omega_g')

        # Neutrino -- Ultrarelativisitc
        self.sym['rho_n_ur'] = sym.symbols('rho_n_ur')
        self.sym['rhohat_n_ur'] = sym.symbols('\hat{rho}_n_ur')
        self.sym['Omega_n_ur'] = sym.symbols('Omega_n_ur')
        
        # Neutrino -- non-relativistic
        self.sym['rho_n_nr'] = sym.symbols('rho_n_nr')
        self.sym['rhohat_n_nr'] = sym.symbols('\hat{rho}_n_nr')
        self.sym['Omega_n_nr'] = sym.symbols('Omega_n_nr')
        self.sym['w_n_nr'] = sym.symbols('w_n_nr')
        
        # Baryon
        self.sym['rho_b'] = sym.symbols('rho_b')
        self.sym['rhohat_b'] = sym.symbols('\hat{rho}_b')
        self.sym['Omega_b'] = sym.symbols('Omega_b')
        
        # Cold dark matter
        self.sym['rho_c'] = sym.symbols('rho_c')
        self.sym['rhohat_c'] = sym.symbols('\hat{rho}_c')
        self.sym['Omega_c'] = sym.symbols('Omega_c')

        # Lambda or Dynamical Dark Energy
        self.sym['rho_l'] = sym.symbols('rho_l')
        self.sym['rhohat_l'] = sym.symbols('\hat{rho}_l')
        self.sym['Omega_l'] = sym.symbols('Omega_l')
        self.sym['w_l'] = sym.symbols('w_l') # Lambda or Dynamical Dark Energy EoS

        self.sym['H0'] = sym.symbols('H_0')

        self.sym['f_phi'] = sym.symbols('f_phi')
        self.sym['Theta'] = sym.symbols('Theta')

        # Symbolic expression dictionary
        self.symfunc = {}

        # Lambda function dictionary
        self.lambda_funcs = {}

        # Parameter values
        self.params['mass_ratios'] = None

        # newton tolerance
        self.newton_tol = None
    
    
    # Check symbolic symbol dictionary

    def _check_sym_key(self, key):
        """
        Checks and returns a boolean to indicate whether a key exists in the dictionary for 
        symbolic variables.

        Parameters
        ----------
        key : str
            Dictionary variable.
        """
        if key in self.sym:
            return True
        else:
            return False
    

    def _check_sym_keys(self, keys):
        """
        Checks and returns a boolean to indicate whether the lists of keys exists in the 
        dictionary for symbolic variables.

        Parameters
        ----------
        keys : list
            List of dictionary variable.
        """
        if all(key in self.sym for key in keys):
            return True
        else:
            return False
    
    # Check symbolic function dictionary

    def _check_symfunc_key(self, key):
        """
        Checks and returns a boolean to indicate whether a key exists in the dictionary for 
        symbolic variables.

        Parameters
        ----------
        key : str
            Dictionary variable.
        """
        if key in self.symfunc:
            return True
        else:
            return False
    

    def _check_symfunc_keys(self, keys):
        """
        Checks and returns a boolean to indicate whether the lists of keys exists in the 
        dictionary for symbolic variables.

        Parameters
        ----------
        keys : list
            List of dictionary variable.
        """
        if all(key in self.symfunc for key in keys):
            return True
        else:
            return False
    
    # Printing functionalities
    
    def get_latex(self, variable, simplify=False):
        """
        Returns the latex string for a given variable function.

        Parameters
        ----------
        variable : str
            Variable for symbolic function computed within the class.
        simplify : bool, optional
            Instructs the code to simplify the expression.
        
        Returns
        -------
        latex_str : str
            Latex string for given variable function.
        """
        if simplify:
            if variable in self.sym:
                variable_sym = self.sym[variable]
                variable_sym_simple = sym.simplify(variable_sym)
                latex_str = sym.latex(variable_sym_simple)
                return latex_str
            elif variable in self.symfunc:
                variable_sym = self.symfunc[variable]
                variable_sym_simple = sym.simplify(variable_sym)
                latex_str = sym.latex(variable_sym_simple)
                return latex_str
            else:
                assert False, 'Symbol/function is currently undefined.'
        else:
            if variable in self.sym:
                variable_sym = self.sym[variable]
                latex_str = sym.latex(variable_sym)
                return latex_str
            elif variable in self.symfunc:
                variable_sym = self.symfunc[variable]
                latex_str = sym.latex(variable_sym)
                return latex_str
            else:
                assert False, 'Symbol/function is currently undefined.'


    def initialise_sympy_print(self):
        """
        Prints functions using pretty-printing.
        """
        sym.init_printing()

    # Horndeski function definitions

    def define_K(self, K_func, K_sym):
        """
        Define the equation for K and symbols.

        Parameters
        ----------
        K_exp : str
            The expression for K in string format.
        K_sym : str
            The symbols in the K function.

        Example
        -------
        For ESS model we would do the following
        self.define_K("K_1*X + K_2*X*X", "K_1 K_2")
        """
        if K_sym is not None:
            _K_sym = sym.symbols(K_sym)
            if isinstance(_K_sym, sym.Symbol):
                self.sym['K_syms'] = [_K_sym]
            elif isinstance(_K_sym, (tuple, list)):
                self.sym['K_syms'] = _K_sym
        else:
            self.sym['K_syms'] = []
        self.symfunc['K'] = sym.sympify(K_func)
    

    def define_G3(self, G3_func, G3_sym):
        """
        Define the equation for K and symbols.

        Parameters
        ----------
        G3_exp : str
            The expression for G3 in string format.
        G3_sym : str
            The symbols in the G3 function.
        """
        if G3_sym is not None:
            _G3_sym = sym.symbols(G3_sym)
            if isinstance(_G3_sym, sym.Symbol):
                self.sym['G3_syms'] = [_G3_sym]
            elif isinstance(_G3_sym, (tuple, list)):
                self.sym['G3_syms'] = _G3_sym
        else:
            self.sym['G3_syms'] = []
        self.symfunc['G3'] = sym.sympify(G3_func)
    
    
    def define_G4(self, G4_func, G4_sym):
        """
        Define the equation for K and symbols.

        Parameters
        ----------
        G4_exp : str
            The expression for G4 in string format.
        G4_sym : str
            The symbols in the G4 function.
        """
        if G4_sym is not None:
            _G4_sym = sym.symbols(G4_sym)
            if isinstance(_G4_sym, sym.Symbol):
                self.sym['G4_syms'] = [_G4_sym]
            elif isinstance(_G4_sym, (tuple, list)):
                self.sym['G4_syms'] = _G4_sym
        else:
            self.sym['G4_syms'] = []
        self.symfunc['G4'] = sym.sympify(G4_func)


    def _get_K_G3_G4_syms(self):
        """
        Combines the symbols for 
        """
        self.sym['K_G3_G4_syms'] = []
        for K_sym in self.sym['K_syms']:
            self.sym['K_G3_G4_syms'].append(K_sym)
        for G3_sym in self.sym['G3_syms']:
            self.sym['K_G3_G4_syms'].append(G3_sym)
        for G4_sym in self.sym['G4_syms']:
            self.sym['K_G3_G4_syms'].append(G4_sym)


    def get_K_derivatives(self):
        """
        Computes the symbolic differentials of the K function with respect to X and phi.
        """
        self.symfunc['Kx'] = sym.diff(self.symfunc['K'],self.sym['X'])
        self.symfunc['Kxx'] = sym.diff(self.symfunc['Kx'],self.sym['X'])
        self.symfunc['Kxphi'] = sym.diff(self.symfunc['Kx'],self.sym['phi'])
        self.symfunc['Kphi'] = sym.diff(self.symfunc['K'],self.sym['phi'])
        self.symfunc['Kphiphi'] = sym.diff(self.symfunc['Kphi'],self.sym['phi'])
        self.symfunc['Kphix'] = sym.diff(self.symfunc['Kphi'],self.sym['X'])
            

    def get_G3_derivatives(self):
        """
        Computes the symbolic differentials of the G3 function with respect to X and phi.
        """
        self.symfunc['G3x'] = sym.diff(self.symfunc['G3'],self.sym['X'])
        self.symfunc['G3xx'] = sym.diff(self.symfunc['G3x'],self.sym['X'])
        self.symfunc['G3xphi'] = sym.diff(self.symfunc['G3x'],self.sym['phi'])
        self.symfunc['G3phi'] = sym.diff(self.symfunc['G3'],self.sym['phi'])
        self.symfunc['G3phiphi'] = sym.diff(self.symfunc['G3phi'],self.sym['phi'])
        self.symfunc['G3phix'] = sym.diff(self.symfunc['G3phi'],self.sym['X'])
            

    def get_G4_derivatives(self):
        """
        Computes the symbolic differentials of the G4 function with respect to X and phi.
        """
        self.symfunc['G4phi'] = sym.diff(self.symfunc['G4'],self.sym['phi'])
        self.symfunc['G4phiphi'] = sym.diff(self.symfunc['G4phi'],self.sym['phi'])

    # Symbolic function definitions
    
    def get_rho_phi(self):
        """
        Returns the scalar field fractional density function, following equation 2.4 in https://arxiv.org/pdf/2209.01666.pdf.
        """
        self.symfunc['rho_phi'] = (1/(2*self.sym['M_G4p2']*self.symfunc['G4']) - 1)
        self.symfunc['rho_phi'] *= (self.sym['rho_g'] + self.sym['rho_n_ur'] + self.sym['rho_n_nr'] + self.sym['rho_b'] + self.sym['rho_c'] + self.sym['rho_l'])/((self.sym['E']**2))
        term1 = self.sym['M_Kp2']*self.sym['X']*self.symfunc['Kx']/(self.sym['E']**2)
        term1 -= self.sym['M_Kp2']*self.symfunc['K']/(2*self.sym['E']**2)
        term1 += 3*self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.sym['phi_prime']*self.symfunc['G3x']
        term1 -= self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.symfunc['G3phi']/(self.sym['E']**2)
        term1 -= 3*self.sym['M_G4p2']*self.sym['phi_prime']*self.symfunc['G4phi']
        self.symfunc['rho_phi'] += (1/(3*self.sym['M_G4p2']*self.symfunc['G4']))*term1
        self.symfunc['rho_phi'] *= self.sym['E']**2

    def get_Omega_phi(self):
        """
        Returns the scalar field fractional density function, following equation 2.4 in https://arxiv.org/pdf/2209.01666.pdf.
        """
        self.symfunc['Omega_phi'] = self.symfunc['rho_phi']/(self.sym['E']**2)


    def get_fried_closure(self):
        """
        Returns the Friedmann closure relation -> equation 2.3 in https://arxiv.org/abs/2209.01666, should be equal to 0. 
        """
        self.get_rho_phi()
        self.symfunc['fried_closure'] = self.sym['rho_g'] + self.sym['rho_n_ur'] + self.sym['rho_n_nr']
        self.symfunc['fried_closure'] += self.sym['rho_b'] + self.sym['rho_c'] + self.sym['rho_l']
        self.symfunc['fried_closure'] += self.symfunc['rho_phi']
        self.symfunc['fried_closure'] -= self.sym['E']**2
        self.symfunc['fried_closure'] /= self.sym['E']**2


    def get_G_G_4(self):
        """
        Returns the G_G_4/G_N function.
        """
        self.symfunc['G_G_4/G_N'] = 1/(2*self.sym['M_G4p2']*self.symfunc['G4'])
    

    def get_A(self):
        """
        Code equivalent of equation 2.7 in https://arxiv.org/abs/2209.01666.
        """
        self.symfunc['A'] = self.sym['M_Kp2']*self.symfunc['Kx']
        self.symfunc['A'] -= 2*self.sym['M_G3p']*self.sym['M_sp']*self.symfunc['G3phi']
        self.symfunc['A'] += 2*self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.symfunc['G3phix']
        term1 = 6*self.sym['M_G3p']*self.sym['M_sp']*self.sym['phi_prime']*(self.symfunc['G3x']+self.sym['X']*self.symfunc['G3xx'])
        term1 += (self.sym['phi_prime']**2)*(self.sym['M_Kp2']*self.symfunc['Kxx']-2*self.sym['M_G3p']*self.sym['M_sp']*self.symfunc['G3phix'])
        term1 *= self.sym['E']**2
        self.symfunc['A'] += term1


    def get_B1(self):
        """
        Code equivalent of equation 2.9 in https://arxiv.org/abs/2209.01666.
        """
        self.symfunc['B1'] = 6*self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.symfunc['G3x'] - 6*self.sym['M_G4p2']*self.symfunc['G4phi']


    def get_B2(self):
        """
        Code equivalent of equation 2.10 in https://arxiv.org/abs/2209.01666.
        """
        term1 = self.sym['M_Kp2']*self.symfunc['Kx']
        term1 -= 2*self.sym['M_G3p']*self.sym['M_sp']*self.symfunc['G3phi']
        term1 += 2*self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.symfunc['G3xphi']
        term1 *= 3*self.sym['phi_prime']
        self.symfunc['B2'] = term1
        term2 = self.sym['M_Kp2']*self.symfunc['Kxphi'] - 2*self.sym['M_G3p']*self.sym['M_sp']*self.symfunc['G3phiphi']
        term2 *= self.sym['phi_prime']**2
        self.symfunc['B2'] += term2
        self.symfunc['B2'] -= self.sym['M_Kp2']*self.symfunc['Kphi']/(self.sym['E']**2)
        self.symfunc['B2'] -= 12*self.sym['M_G4p2']*self.symfunc['G4phi']
        self.symfunc['B2'] += 18*self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.symfunc['G3x']
        self.symfunc['B2'] += 2*self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.symfunc['G3phiphi']/(self.sym['E']**2)


    def get_E_prime(self):
        """
        Returns the evolution of the dimensionless Hubble function E=H/H0, as a function of log(a), follows equation 2.5 
        https://arxiv.org/pdf/2209.01666.pdf modified to remove the dependence on phi_primeprime.
        """
        term1 = - self.sym['M_Kp2']*self.symfunc['K']/(2*self.sym['E']**2)
        term1 += self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.symfunc['G3phi']/(self.sym['E']**2)
        term1 -= 2*self.sym['M_G4p2']*self.sym['phi_prime']*self.symfunc['G4phi']
        term1 -= 2*self.sym['M_G4p2']*self.sym['X']*self.symfunc['G4phiphi']/(2*self.sym['E']**2)
        term1 -= (3/2)*((1/3)*self.sym['rho_g'] + (1/3)*self.sym['rho_n_ur'] + self.sym['w_n_nr']*self.sym['rho_n_nr'] + self.sym['w_l']*self.sym['rho_l'])/(self.sym['E']**2)
        term1 -= 3*self.sym['M_G4p2']*self.symfunc['G4']
        term1 *= self.symfunc['A']
        term2 = self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.symfunc['G3x'] - self.sym['M_G4p2']*self.symfunc['G4phi']
        term2 *= self.symfunc['B2']
        calA = term1 - term2
        calB = 2*self.sym['M_G4p2']*self.symfunc['A']*self.symfunc['G4']
        calB += self.symfunc['B1']*(self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.symfunc['G3x'] - self.sym['M_G4p2']*self.symfunc['G4phi'])
        self.symfunc['E_prime'] = self.sym['E']*calA/calB
        
    
    def get_phi_primeprime(self):
        """
        Returns equation 2.6 to compute phi_primeprime.
        """
        self.symfunc['phi_primeprime'] = - (self.sym['phi_prime'] + self.symfunc['B1']/self.symfunc['A'])*self.symfunc['E_prime']/self.sym['E'] 
        self.symfunc['phi_primeprime'] -= self.symfunc['B2']/self.symfunc['A']


    def get_calE(self):
        """
        Code equivalent of equation 5 in https://arxiv.org/abs/1111.6749, given in units of H0^2 * Mp^2.
        """
        self.symfunc['calE'] = 2*self.sym['M_Kp2']*self.sym['X']*self.symfunc['Kx']
        self.symfunc['calE'] -= self.sym['M_Kp2']*self.symfunc['K']
        self.symfunc['calE'] += 6*self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*(self.sym['E']**2)*self.sym['phi_prime']*self.symfunc['G3x']
        self.symfunc['calE'] -= 2*self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.symfunc['G3phi']
        self.symfunc['calE'] -= 6*self.sym['M_G4p2']*(self.sym['E']**2)*self.symfunc['G4']
        self.symfunc['calE'] -= 6*self.sym['M_G4p2']*(self.sym['E']**2)*self.sym['phi_prime']*self.symfunc['G4phi']


    def get_calP(self):
        """
        Code equivalent of equation 6 in https://arxiv.org/abs/1111.6749, given in units of H0^2 * Mp^2.
        """
        self.symfunc['calP'] = self.sym['M_Kp2']*self.symfunc['K']
        term2_bracket = self.symfunc['G3phi'] + ((self.symfunc['E_prime']/self.sym['E'])*self.sym['phi_prime'] + self.symfunc['phi_primeprime'])*(self.sym['E']**2)*self.symfunc['G3x']
        self.symfunc['calP'] -= 2*self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*term2_bracket
        self.symfunc['calP'] += 2*self.sym['M_G4p2']*(3 + 2*self.symfunc['E_prime']/self.sym['E'])*(self.sym['E']**2)*self.symfunc['G4']
        term4_bracket = ((self.symfunc['E_prime']/self.sym['E'])*self.sym['phi_prime'] + self.symfunc['phi_primeprime'] + 2*self.sym['phi_prime'])
        self.symfunc['calP'] += 2*self.sym['M_G4p2']*term4_bracket*(self.sym['E']**2)*self.symfunc['G4phi']
        self.symfunc['calP'] += 4*self.sym['M_G4p2']*self.sym['X']*self.symfunc['G4phiphi']


    def get_theta(self):
        """
        Code equivalent of equation 2.15 in https://arxiv.org/abs/2209.01666, given in units of H0*Mp^2.
        """
        self.symfunc['theta'] = -self.sym['M_sp']*self.sym['M_G3p']*self.sym['phi_prime']*self.sym['X']*self.symfunc['G3x']
        self.symfunc['theta'] += 2*self.sym['M_G4p2']*self.symfunc['G4']
        self.symfunc['theta'] += self.sym['M_G4p2']*self.sym['phi_prime']*self.symfunc['G4phi']
        self.symfunc['theta'] *= self.sym['E']

    
    def get_alpha0(self):
        """
        Code equivalent of equation 3.5 (see 2.16-2.19) in https://arxiv.org/abs/2209.01666.
        Given by alpha0 = A0/2G4.
        """

        self.get_calE() # units H0^{2} * Mp^2
        self.get_calP() # units H0^{2} * Mp^2
        self.get_theta() # units H0 * Mp^2

        X_prime = (self.sym['E']**2)*self.sym['phi_prime']*(self.symfunc['phi_primeprime'] + self.symfunc['E_prime']*self.sym['phi_prime']/self.sym['E'])

        theta_prime = self.symfunc['E_prime']*sym.diff(self.symfunc['theta'], self.sym['E'])
        theta_prime += self.sym['phi_prime']*sym.diff(self.symfunc['theta'], self.sym['phi'])
        theta_prime += X_prime*sym.diff(self.symfunc['theta'], self.sym['X'])

        self.symfunc['alpha0'] = theta_prime/self.sym['E']
        self.symfunc['alpha0'] += self.symfunc['theta']/self.sym['E']
        self.symfunc['alpha0'] -= 2*self.sym['M_G4p2']*self.symfunc['G4']
        self.symfunc['alpha0'] -= 4*self.sym['M_G4p2']*self.sym['phi_prime']*self.symfunc['G4phi']
        self.symfunc['alpha0'] -= (self.symfunc['calE'] + self.symfunc['calP'])/(2*self.sym['E']**2)
        self.symfunc['alpha0'] /= 2*self.sym['M_G4p2']*self.symfunc['G4']
        

    def get_alpha1(self):
        """
        Code equivalent of equation 3.5 (see 2.16-2.19) in https://arxiv.org/abs/2209.01666.
        Given by alpha1 = A1/2G4.
        """
        self.symfunc['alpha1'] = self.sym['phi_prime']*self.symfunc['G4phi']/self.symfunc['G4']


    def get_alpha2(self):
        """
        Code equivalent of equation 3.5 (see 2.16-2.19) in https://arxiv.org/abs/2209.01666.
        Given by alpha2 = A2/2G4.
        """
        self.get_theta() # units H0 * Mp^2
        self.symfunc['alpha2'] = 2*self.sym['M_G4p2']*self.symfunc['G4'] - self.symfunc['theta']/self.sym['E']
        self.symfunc['alpha2'] /= 2*self.sym['M_G4p2']*self.symfunc['G4']


    def get_beta0(self):
        """
        Code equivalent of equation 3.5 (see 2.16-2.19) in https://arxiv.org/abs/2209.01666.
        Given by beta0 = B0/2G4.
        """
        self.symfunc['beta0'] = self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.sym['phi_prime']*self.symfunc['G3x']
        self.symfunc['beta0'] /= 2*self.sym['M_G4p2']*self.symfunc['G4']
    

    def get_calB(self):
        """
        Code equivalent of equation 3.7 in https://arxiv.org/abs/2209.01666.
        """
        self.get_alpha0()
        self.get_alpha1()
        self.get_alpha2()
        self.get_beta0()
        self.symfunc['calB'] = 4.*self.symfunc['beta0']/(self.symfunc['alpha0'] + 2.*self.symfunc['alpha1']*self.symfunc['alpha2'] + self.symfunc['alpha2']*self.symfunc['alpha2'])


    def get_calC(self):
        """
        Code equivalent of equation 3.7 in https://arxiv.org/abs/2209.01666.
        """
        self.get_alpha0()
        self.get_alpha1()
        self.get_alpha2()
        self.symfunc['calC'] = (self.symfunc['alpha1'] + self.symfunc['alpha2'])/(self.symfunc['alpha0'] + 2.*self.symfunc['alpha1']*self.symfunc['alpha2'] + self.symfunc['alpha2']*self.symfunc['alpha2'])


    def get_beta(self):
        """
        The coupling in equation 3.13 in https://arxiv.org/abs/2209.01666, i.e. the deviation from 1.
        """
        self.get_alpha1()
        self.get_alpha2()
        self.get_calC()
        self.symfunc['beta'] = -1.*(self.symfunc['alpha1'] + self.symfunc['alpha2'])*self.symfunc['calC']
    
    
    # # Symbolic definitions for Bellini & Sawicki alphas

    def get_M_star_sq(self):
        """
        Code equivalent of equation A.6 in https://iopscience.iop.org/article/10.1088/1475-7516/2014/07/050 within reduced Horndeski class
        given in units of Mp^2
        """
        self.symfunc['M_star_sq'] = 2*self.sym['M_G4p2']*self.symfunc['G4']
    

    def get_alpha_M(self):
        """
        Code equivalent of equation A.7 in https://iopscience.iop.org/article/10.1088/1475-7516/2014/07/050 within reduced Horndeski class.
        """
        self.symfunc['alpha_M'] = 2*self.sym['M_G4p2']*self.sym['phi_prime']*self.symfunc['G4phi'] # in units of Mp^2
        self.symfunc['alpha_M'] /= self.symfunc['M_star_sq']
    

    def get_alpha_K(self):
        """
        Code equivalent of equation A.8 in https://iopscience.iop.org/article/10.1088/1475-7516/2014/07/050 within reduced Horndeski class.
        """
        term1 = self.sym['M_Kp2']*self.symfunc['Kx'] 
        term1 += 2*self.sym['M_Kp2']*self.sym['X']*self.symfunc['Kxx']
        term1 -= 2*self.sym['M_G3p']*self.sym['M_sp']*self.symfunc['G3phi']
        term1 -= 2*self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.symfunc['G3phix']
        self.symfunc['alpha_K'] = 2*self.sym['X']*term1
        self.symfunc['alpha_K'] += 12*(self.sym['E']**2)*self.sym['M_G3p']*self.sym['M_sp']*self.sym['phi_prime']*self.sym['X']*(self.symfunc['G3x']+self.sym['X']*self.symfunc['G3xx'])
        self.symfunc['alpha_K'] /= self.sym['E']**2 # in units of Mp^2
        self.symfunc['alpha_K'] /= self.symfunc['M_star_sq']


    def get_alpha_B(self):
        """
        Code equivalent of equation A.9 in https://iopscience.iop.org/article/10.1088/1475-7516/2014/07/050 within reduced Horndeski class.
        """
        self.symfunc['alpha_B'] = 2*self.sym['phi_prime']
        self.symfunc['alpha_B'] *= self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.symfunc['G3x'] - self.sym['M_G4p2']*self.symfunc['G4phi'] # in units of M_p^{2}
        self.symfunc['alpha_B'] /= self.symfunc['M_star_sq']

    # # Symbolic definitions for scalar field density and pressure

    def get_tilde_calE(self):
        """
        Code equivalent of equation A.1 in https://iopscience.iop.org/article/10.1088/1475-7516/2014/07/050 within reduced Horndeski class.
        Referred to as Tilde epsilon in internal Hi-COLA notes. Given in units of H0^2.
        """
        self.symfunc['tilde_calE'] = -self.sym['M_Kp2']*self.symfunc['K']
        self.symfunc['tilde_calE'] += 2*self.sym['M_Kp2']*self.sym['X']*self.symfunc['Kx']
        self.symfunc['tilde_calE'] -= 2*self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.symfunc['G3phi']
        term1 = self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.symfunc['G3x'] - self.sym['M_G4p2']*self.symfunc['G4phi']
        self.symfunc['tilde_calE'] += 6*(self.sym['E']**2)*self.sym['phi_prime']*term1
        self.symfunc['tilde_calE'] /= self.symfunc['M_star_sq']
    

    def get_tilde_calP(self):
        """
        Code equivalent of equation A.2 in https://iopscence.iop.org/article/10.1088/1475-7516/2014/07/050 within reduced Horndeski class.
        Referred to as Tilde P in internal Hi-COLA notes. Given in units of H0^2.
        """
        self.symfunc['tilde_calP'] = self.sym['M_Kp2']*self.symfunc['K']
        self.symfunc['tilde_calP'] -= 2*self.sym['X']*(self.sym['M_G3p']*self.sym['M_sp']*self.symfunc['G3phi'] - 2*self.sym['M_G4p2']*self.symfunc['G4phiphi'])
        self.symfunc['tilde_calP'] += 4*self.sym['M_G4p2']*(self.sym['E']**2)*self.sym['phi_prime']*self.symfunc['G4phi']
        term1 = self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.symfunc['G3x'] - self.sym['M_G4p2']*self.symfunc['G4phi']
        term1 *= self.symfunc['E_prime']*self.sym['phi_prime']/self.sym['E'] + self.symfunc['phi_primeprime']
        term1 *= 2*(self.sym['E']**2)
        self.symfunc['tilde_calP'] -= term1
        self.symfunc['tilde_calP'] /= self.symfunc['M_star_sq']


    # # Stability and sound speed related functions

    def get_Q_s(self):
        """
        Code equivalent of the LHS of the first inequality in 3.13 in https://iopscience.iop.org/article/10.1088/1475-7516/2014/07/050 within reduced Horndeski class.
        Q_s is given in units of Mp^2.
        """
        self.symfunc['D'] = self.symfunc['alpha_K'] + (3/2)*(self.symfunc['alpha_B']**2)
        self.symfunc['Q_s'] = 2*self.symfunc['M_star_sq']*self.symfunc['D']/((2-self.symfunc['alpha_B'])**2) 
    

    def get_c_s_sq(self):
        """
        Code equivalent of the LHS of the second inequality in 3.13 in https://iopscience.iop.org/article/10.1088/1475-7516/2014/07/050 within reduced Horndeski class.
        """
        alpha_B = self.symfunc['alpha_B'].subs({self.sym['X']: (self.sym['E']**2 * self.sym['phi_prime']**2)/2})

        alpha_B_prime = self.symfunc['E_prime']*sym.diff(alpha_B, self.sym['E'])
        alpha_B_prime += self.sym['phi_prime']*sym.diff(alpha_B, self.sym['phi'])
        alpha_B_prime += self.symfunc['phi_primeprime']*sym.diff(alpha_B, self.sym['phi_prime'])
        self.symfunc['alpha_B_prime'] = alpha_B_prime

        self.symfunc['c_s_sq_D'] = 2-self.symfunc['alpha_B']
        self.symfunc['c_s_sq_D'] *= (self.symfunc['E_prime']/self.sym['E'] - 0.5*self.symfunc['alpha_B'] - self.symfunc['alpha_M'])
        self.symfunc['c_s_sq_D'] -= alpha_B_prime
        term1 = self.sym['rho_g']*(1 + 1/3)
        term1 += self.sym['rho_n_ur']*(1 + 1/3)
        term1 += self.sym['rho_n_nr']*(1 + self.sym['w_n_nr'])
        term1 += self.sym['rho_b']
        term1 += self.sym['rho_c']
        term1 += self.sym['rho_l']*(1 + self.sym['w_l'])
        self.symfunc['c_s_sq_D'] += 3*(term1/(self.sym['E']**2))/self.symfunc['M_star_sq']
        self.symfunc['c_s_sq_D'] *= -1
        self.symfunc['c_s_sq'] = self.symfunc['c_s_sq_D']/self.symfunc['D']
    

    # Set mass ratios
    
    def set_mass_ratios(self, M_sp=1, M_Kp2=1, M_G3p=1, M_G4p2=1):
        """
        Allows the users to assign specific values to mass ratios, those unchanged will be set to 1.

        Parameters
        ----------
        M_sp : float, optional
            Mass associated with the scalar field with respect to the Planck mass.
        M_Kp2 : float, optional
            Mass associated with the Kinetic expression with respect to the Planck mass squared.
        M_G3p : float, optional
            Mass associated with the G3 expression with respect to the Planck mass squared.
        M_G4p : float, optional
            Mass associated with the G4 expression with respect to the Planck mass squared.
        """
        self.params['mass_ratios'] = {
            'M_sp': M_sp, 'M_Kp2': M_Kp2, 'M_G3p': M_G3p, 'M_G4p2': M_G4p2
        }
    

    def get_absolute_mass_ratios(self):
        """
        Allows the users to assign specific values to mass ratios, those unchanged will be set to 1.

        Parameters
        ----------
        M_p : float, optional
            Planck mass.
        M_pG4 : float, optional
            Mass ratio M_p/M_G4.
        M_KG4 : float, optional
            Mass ratio M_K/M_G4.
        M_G3s : float, optional
            Mass ratio M_G3/M_s.
        M_sG4 : float, optional
            Mass ratio M_s/M_G4.
        M_G3G4 : float, optional
            Mass ratio M_G4/M_G4.
        M_Ks : float, optional
            Mass ratio M_K/M_s.
        M_gp : float, optional
            Mass ratio M_g/M_p.
        """
        if self.params['mass_ratios'] is None:
            self.set_mass_ratios()
    

    # Construct the symbolic model

    def construct_model(self):
        """
        Constructs Horndeski model with user defined functions.
        """

        if self.verbose:
            print('Hi-COLA: Constructing model')

        if self._check_symfunc_keys(['K', 'G3', 'G4']) == False:
            assert False, 'Functions for K, G3 and G4 remain undefined.'
        else:
            self._get_K_G3_G4_syms()
            self.get_K_derivatives()
            self.get_G3_derivatives()
            self.get_G4_derivatives()

            self.get_absolute_mass_ratios()

            if self.verbose:
                print(' - substituting Ehat = E*fH, phihat = phi/fH, phihat_prime = phi_prime/fH, rhohat_i = rho_i*fH**2')
                
            sub_dict = {
                self.sym['E']: self.sym['E']*self.sym['f_H'],
                self.sym['phi']: self.sym['phi']/self.sym['f_H'],
                self.sym['phi_prime']: self.sym['phi_prime']/self.sym['f_H'],
                self.sym['rho_g']: self.sym['rho_g']*self.sym['f_H']**2,
                self.sym['rho_b']: self.sym['rho_b']*self.sym['f_H']**2,
                self.sym['rho_c']: self.sym['rho_c']*self.sym['f_H']**2,
                self.sym['rho_l']: self.sym['rho_l']*self.sym['f_H']**2,
                self.sym['rho_n_ur']: self.sym['rho_n_ur']*self.sym['f_H']**2,
                self.sym['rho_n_nr']: self.sym['rho_n_nr']*self.sym['f_H']**2,
            }

            if self.verbose:
                print(' - into symbolic functions...')

            self.get_G_G_4()
            G_G_4_GN = self.symfunc['G_G_4/G_N'].subs(sub_dict)

            self.get_rho_phi()
            rho_phi = self.symfunc['rho_phi'].subs(sub_dict)

            self.get_Omega_phi()
            Omega_phi = self.symfunc['Omega_phi'].subs(sub_dict)

            self.get_fried_closure()
            fried_closure = self.symfunc['fried_closure'].subs(sub_dict)

            self.get_A()
            A = self.symfunc['A'].subs(sub_dict)

            self.get_B1()
            B1 = self.symfunc['B1'].subs(sub_dict)

            self.get_B2()
            B2 = self.symfunc['B2'].subs(sub_dict)

            self.get_E_prime()
            E_prime = self.symfunc['E_prime'].subs(sub_dict)

            self.get_phi_primeprime()
            phi_primeprime = self.symfunc['phi_primeprime'].subs(sub_dict)

            self.get_alpha0()
            alpha0 = self.symfunc['alpha0'].subs(sub_dict)

            self.get_alpha1()
            alpha1 = self.symfunc['alpha1'].subs(sub_dict)

            self.get_alpha2()
            alpha2 = self.symfunc['alpha2'].subs(sub_dict)

            self.get_beta0()
            beta0 = self.symfunc['beta0'].subs(sub_dict)

            self.get_calB()
            calB = self.symfunc['calB'].subs(sub_dict)

            self.get_calC()
            calC = self.symfunc['calC'].subs(sub_dict)

            self.get_beta()
            beta = self.symfunc['beta'].subs(sub_dict)

            self.get_M_star_sq()
            M_star_sq = self.symfunc['M_star_sq'].subs(sub_dict)
            
            self.get_alpha_M()
            alpha_M = self.symfunc['alpha_M'].subs(sub_dict)

            self.get_alpha_B()
            alpha_B = self.symfunc['alpha_B'].subs(sub_dict)

            self.get_alpha_K()
            alpha_K = self.symfunc['alpha_K'].subs(sub_dict)
            
            self.get_tilde_calE()
            tilde_calE = self.symfunc['tilde_calE'].subs(sub_dict)

            self.get_tilde_calP()
            tilde_calP = self.symfunc['tilde_calP'].subs(sub_dict)

            self.get_Q_s()
            D = self.symfunc['D'].subs(sub_dict)
            Q_s = self.symfunc['Q_s'].subs(sub_dict)

            self.get_c_s_sq()
            alpha_B_prime = self.symfunc['alpha_B_prime'].subs(sub_dict)
            c_s_sq_D = self.symfunc['c_s_sq_D'].subs(sub_dict)
            c_s_sq = self.symfunc['c_s_sq'].subs(sub_dict)

            if self.verbose:
                print(' - substituting X = 0.5 * E^2 * phi_prime^2')
                print(
                    ' - substituting mass ratios: M_sp=%0.2f, M_Kp2=%0.2f, M_G3p=%0.2f, M_G4p2=%0.2f' % (
                        self.params['mass_ratios']['M_sp'], self.params['mass_ratios']['M_Kp2'], 
                        self.params['mass_ratios']['M_G3p'], self.params['mass_ratios']['M_G4p2'],
                    )
                )
            
            Xreal = 0.5*(self.sym['E']**2)*(self.sym['phi_prime']**2.)

            sub_dict = {
                self.sym['X']: Xreal,
                self.sym['M_sp']: self.params['mass_ratios']['M_sp'],
                self.sym['M_Kp2']: self.params['mass_ratios']['M_Kp2'],
                self.sym['M_G3p']: self.params['mass_ratios']['M_G3p'],
                self.sym['M_G4p2']: self.params['mass_ratios']['M_G4p2'],
            }

            if self.verbose:
                print(' - into symbolic functions...')

            G_G_4_GN = G_G_4_GN.subs(sub_dict)

            rho_phi = rho_phi.subs(sub_dict)

            Omega_phi = Omega_phi.subs(sub_dict)

            fried_closure = fried_closure.subs(sub_dict)

            fried_closure_dE = sym.diff(fried_closure, self.sym['E']).subs(sub_dict)

            fried_closure_dE2 = sym.diff(fried_closure_dE, self.sym['E']).subs(sub_dict)

            A = A.subs(sub_dict)

            B1 = B1.subs(sub_dict)

            B2 = B2.subs(sub_dict)

            E_prime = E_prime.subs(sub_dict)

            phi_primeprime = phi_primeprime.subs(sub_dict)

            alpha0 = alpha0.subs(sub_dict)

            alpha1 = alpha1.subs(sub_dict)

            alpha2 = alpha2.subs(sub_dict)

            beta0 = beta0.subs(sub_dict)

            calB = calB.subs(sub_dict)

            calC = calC.subs(sub_dict)

            beta = beta.subs(sub_dict)

            M_star_sq = M_star_sq.subs(sub_dict)
            
            alpha_M = alpha_M.subs(sub_dict)

            alpha_B = alpha_B.subs(sub_dict)

            alpha_K = alpha_K.subs(sub_dict)
            
            tilde_calE = tilde_calE.subs(sub_dict)

            tilde_calP = tilde_calP.subs(sub_dict)

            D = D.subs(sub_dict)
            Q_s = Q_s.subs(sub_dict)

            c_s_sq_D = c_s_sq_D.subs(sub_dict)
            alpha_B_prime = alpha_B_prime.subs(sub_dict)
            c_s_sq = c_s_sq.subs(sub_dict)

            # we will copy these substituted and simplified functions to the class symfunc dictionary, we avoided 
            # doing this before as some of these are re-called and redefined in the 'get' functions.
            self.symfunc['G_G_4/G_N'] = G_G_4_GN
            self.symfunc['rho_phi'] = rho_phi
            self.symfunc['Omega_phi'] = Omega_phi

            self.symfunc['fried_closure'] = fried_closure
            self.symfunc['fried_closure_dE'] = fried_closure_dE
            self.symfunc['fried_closure_dE2'] = fried_closure_dE2
            self.symfunc['A'] = A
            self.symfunc['B1'] = B1
            self.symfunc['B2'] = B2
            self.symfunc['E_prime'] = E_prime
            self.symfunc['phi_primeprime'] = phi_primeprime

            self.symfunc['alpha0'] = alpha0
            self.symfunc['alpha1'] = alpha1
            self.symfunc['alpha2'] = alpha2
            self.symfunc['beta0'] = beta0
            self.symfunc['calB'] = calB
            self.symfunc['calC'] = calC
            self.symfunc['beta'] = beta

            self.symfunc['M_star_sq'] = M_star_sq
            self.symfunc['alpha_M'] = alpha_M
            self.symfunc['alpha_B'] = alpha_B
            self.symfunc['alpha_B_prime'] = alpha_B_prime
            self.symfunc['alpha_K'] = alpha_K
            self.symfunc['tilde_calE'] = tilde_calE
            self.symfunc['tilde_calP'] = tilde_calP
            self.symfunc['D'] = D
            self.symfunc['Q_s'] = Q_s
            self.symfunc['c_s_sq_D'] = c_s_sq_D
            self.symfunc['c_s_sq'] = c_s_sq

            # Lambdify functions

            if self.verbose:
                print(" - 'Lambdify'ing symbolic functions")

            # keep variables fixed to functions solved in the ODE + Horndeski variables
            variables = [
                self.sym['E'],
                self.sym['phi'],
                self.sym['phi_prime'],
                self.sym['rho_g'],
                self.sym['rho_b'],
                self.sym['rho_c'],
                self.sym['rho_l'],
                self.sym['rho_n_ur'],
                self.sym['rho_n_nr'],
                self.sym['w_n_nr'],
                self.sym['w_l'],
                *self.sym['K_G3_G4_syms'],
                self.sym['f_H']
            ]

            # G_G_4/G_N function
            self.lambda_funcs['G_G_4/G_N'] = sym.lambdify(variables, self.symfunc['G_G_4/G_N'], 'numpy')

            # The effective scalar field density
            self.lambda_funcs['rho_phi'] = sym.lambdify(variables, self.symfunc['rho_phi'], 'numpy')
            self.lambda_funcs['Omega_phi'] = sym.lambdify(variables, self.symfunc['Omega_phi'], 'numpy')

            # Friedmann closure relation
            self.lambda_funcs['fried_closure'] = sym.lambdify(variables, self.symfunc['fried_closure'], 'numpy')

            # Friedmann closure relation differentiated by E
            self.lambda_funcs['fried_closure_dE'] = sym.lambdify(variables, self.symfunc['fried_closure_dE'], 'numpy')
            self.lambda_funcs['fried_closure_dE2'] = sym.lambdify(variables, self.symfunc['fried_closure_dE2'], 'numpy')

            # E_prime functions
            self.lambda_funcs['A'] = sym.lambdify(variables, self.symfunc['A'], 'numpy')
            self.lambda_funcs['B1'] = sym.lambdify(variables, self.symfunc['B1'], 'numpy')
            self.lambda_funcs['B2'] = sym.lambdify(variables, self.symfunc['B2'], 'numpy')
            self.lambda_funcs['E_prime'] = sym.lambdify(variables, self.symfunc['E_prime'], 'numpy')

            # phi_primeprime functions
            self.lambda_funcs['phi_primeprime'] = sym.lambdify(variables, self.symfunc['phi_primeprime'], 'numpy')

            # internal Hi-COLA functions for force coupling
            self.lambda_funcs['alpha0'] = sym.lambdify(variables, self.symfunc['alpha0'], 'numpy')
            self.lambda_funcs['alpha1'] = sym.lambdify(variables, self.symfunc['alpha1'], 'numpy')
            self.lambda_funcs['alpha2'] = sym.lambdify(variables, self.symfunc['alpha2'], 'numpy')
            self.lambda_funcs['beta0'] = sym.lambdify(variables, self.symfunc['beta0'], 'numpy')
            self.lambda_funcs['calB'] = sym.lambdify(variables, self.symfunc['calB'], 'numpy')
            self.lambda_funcs['calC'] = sym.lambdify(variables, self.symfunc['calC'], 'numpy')
            self.lambda_funcs['beta'] = sym.lambdify(variables, self.symfunc['beta'], 'numpy')
            
            # Bellini alphas and stability equations
            self.lambda_funcs['M_star_sq'] = sym.lambdify(variables, self.symfunc['M_star_sq'], 'numpy')
            self.lambda_funcs['alpha_M'] = sym.lambdify(variables, self.symfunc['alpha_M'], 'numpy')
            self.lambda_funcs['alpha_B'] = sym.lambdify(variables, self.symfunc['alpha_B'], 'numpy')
            self.lambda_funcs['alpha_B_prime'] = sym.lambdify(variables, self.symfunc['alpha_B_prime'], 'numpy')
            self.lambda_funcs['alpha_K'] = sym.lambdify(variables, self.symfunc['alpha_K'], 'numpy')
            self.lambda_funcs['tilde_calE'] = sym.lambdify(variables, self.symfunc['tilde_calE'], 'numpy')
            self.lambda_funcs['tilde_calP'] = sym.lambdify(variables, self.symfunc['tilde_calP'], 'numpy')
            self.lambda_funcs['D'] = sym.lambdify(variables, self.symfunc['D'], 'numpy')
            self.lambda_funcs['Q_s'] = sym.lambdify(variables, self.symfunc['Q_s'], 'numpy')
            self.lambda_funcs['c_s_sq_D'] = sym.lambdify(variables, self.symfunc['c_s_sq_D'], 'numpy')
            self.lambda_funcs['c_s_sq'] = sym.lambdify(variables, self.symfunc['c_s_sq'], 'numpy')

            if self.verbose:
                print(' - Done!')
    

    def set_cosmo_params(self, H0_ref, Omega_c0_ref, Omega_b0_ref, fphi, K_G3_G4_values, w0=-1., wa=0., Tcmb=2.7255, Tnu0=1.9518, mnu=[0.,0.,0.], Neff=3.044):
        """
        Set cosmological and Horndeski parameters.

        Parameters
        ----------
        H0_ref : float
            Reference LCDM Hubble constant.
        Omega_c0_ref : float
            Reference LCDM cold dark matter density at redshift zero.
        Omega_b0_ref : float
            Reference LCDM baryon density at redshift zero.
        fphi : float
            The fraction of the scalar field density as a fraction of the full dark energy density (including a cosmological constant).
        K_G3_G4_values : list
            A list of values for the Horndeski specific variables. This must match the length of the user defined variable. 
            Check self.sym['K_G3_G4_syms'] to see what variables are expected.
        w0 : float, optional
            Set dark energy equation of state today.
        wa : float, optional
            Set dark energy equation of state gradient.
        Tcmb : float, optional
            The cmb temperature today, set to 2.7255.
        Tnu0 : float, optional
            The relic neutrino temperature today, set to 1.9518.
        mnu : list, optional
            Neutrino mass for each massive species.
        Neff : float, optional
            Effective number of neutrino species used for computing the drag epoch.
        """
        self.params['H0_ref'] = H0_ref
        self.params['Omega_c0_ref'] = Omega_c0_ref
        self.params['Omega_b0_ref'] = Omega_b0_ref
        # photon density set by CMB temperature
        self.params['Tcmb0'] = Tcmb
        self.params['Tnu0'] = Tnu0
        self._get_C_gamma()
        self.params['Omega_g0_ref'] = self.get_Omega_g(1., 1e-2*self.params['H0_ref'], 1.)
        # neutrino density set by CMB temperature
        self._get_C_nu()
        self.params['mnu'] = mnu
        self.params['N_ur'] = 3-len(mnu)
        self.params['Neff'] = Neff
        self.params['Omega_nu_ur0_ref'] = self.get_Omega_nu_ur(1., 1e-2*self.params['H0_ref'], 1.)
        self.params['Omega_nu_nr0_ref'] = 0.
        for _mnu in self.params['mnu']: 
            self.params['Omega_nu_nr0_ref'] = self.get_Omega_nu_nr(1., 1e-2*self.params['H0_ref'], _mnu, 1.) 
        self.params['fphi'] = fphi
        self.params['Omega_l0_LCDM'] = 1. - self.params['Omega_g0_ref'] - self.params['Omega_nu_ur0_ref'] - self.params['Omega_nu_nr0_ref'] - self.params['Omega_c0_ref'] - self.params['Omega_b0_ref']
        self.params['Omega_phi0_ref'] = fphi*self.params['Omega_l0_LCDM']
        self.params['Omega_l0_ref'] = self.params['Omega_l0_LCDM'] - self.params['Omega_phi0_ref']
        assert len(K_G3_G4_values) == len(self.sym['K_G3_G4_syms']), "Length of Horndeski K_G3_G4_values must match number of defined K, G3, G4 variables."
        self.params['K_G3_G4_values'] = K_G3_G4_values
        self.params['fH'] = 1.
        self.params['w0'] = w0
        self.params['wa'] = wa

    # Initiates solver status

    def _initiate_solver_status(self):
        """
        Initiate solver status monitor.
        """
        self._solver_success = True


    # Computes derivaties for the solver via scipy's solve_ivp method

    def _compute_primes(self, x, Y, timeout=5):
        """
        Compute prime functions for numerical solver.

        Parameters
        ----------
        x : float
            Current value of log(a).
        Y : list
            List containing current [E, phi, phi_prime] values.
        timeout : float, optional
            Time in seconds to force the solver to fail.
        """
        # convert x = log(a) to scale factor
        a = np.exp(x)

        # `_` used to denote current value.
        _phi, _phi_prime = Y

        _w_l = self.compute_w_l(a)
        
        _rho_g = self.get_rho_g(a, self.params['H0_ref']*1e-2)
        _rho_b = self.get_rho_b(a, self.params['Omega_b0_ref'])
        _rho_c = self.get_rho_c(a, self.params['Omega_c0_ref'])
        _rho_l = self.get_rho_l(a, self.params['Omega_l0_ref'], w0=self.params['w0'], wa=self.params['wa'])
        
        _rho_nu_ur = self.get_rho_nu_ur(a, self.params['H0_ref']*1e-2)
        _rho_nu_nr = 0.
        for _mnu in self.params['mnu']:
            _rho_nu_nr += self.get_rho_nu_nr(a, self.params['H0_ref']*1e-2, _mnu)
        _w_nu_nr = 0.
        for _mnu in self.params['mnu']:
            _w_nu_nr += self.compute_w_nu_nr(a, _mnu)*self.get_rho_nu_nr(a, self.params['H0_ref']*1e-2, _mnu)
        if _rho_nu_nr != 0.:
            _w_nu_nr /= _rho_nu_nr

        variables = [
            _phi, _phi_prime, 
            _rho_g, _rho_b, _rho_c, _rho_l, 
            _rho_nu_ur, _rho_nu_nr, _w_nu_nr, _w_l, 
            *self.params['K_G3_G4_values'], self.params['fH']
        ]

        _E_guess = self.compute_E_LCDM(a)

        _E = newton(
            lambda _E: self.lambda_funcs['fried_closure'](_E, *variables), 
            _E_guess,
            fprime = lambda _E: self.lambda_funcs['fried_closure_dE'](_E, *variables), 
            fprime2 = lambda _E: self.lambda_funcs['fried_closure_dE2'](_E, *variables), 
            tol=self.newton_tol
        )

        phi_prime = _phi_prime
        phi_primeprime = self.lambda_funcs['phi_primeprime'](_E, *variables)

        timenow = self._check_timer()

        if timenow >= timeout:
            phi_prime = np.nan
            phi_primeprime = np.nan

        if np.isfinite([phi_prime, phi_primeprime]).all() == False:
            self._solver_success = False
        
        return [phi_prime, phi_primeprime]


    # Numerically computed quantities

    def compute_chi_over_delta(self, a, E, calB, calC, G_G_4_G_N):
        """
        Computes the chi/delta function numerically.

        Parameters
        ----------
        a : float or array
            Scale factor.
        E : float or array
            Normalised Hubble expansion.
        calB : float or array
            Code equivalent of equation 3.7 in https://arxiv.org/abs/2209.01666.
        calC : float or array
            Code equivalent of equation 3.7 in https://arxiv.org/abs/2209.01666.
        G_G_4_G_N : float or array
            G_G_4/G_N is the effective Newton's constant, see equation 3.11 and 3.12 in 
            https://arxiv.org/abs/2209.01666 and explanation provided.

        Return
        ------
        chioverdelta : float or array   
            Code equivalent of equation 3.14 in https://arxiv.org/abs/2209.01666.
        """
        if self.params['Omega_c0'] is None:
            chioverdelta = None
        else:
            chioverdelta = calB * calC * (self.params['Omega_c0']+self.params['Omega_b0'])/((E**2)*(a**3)) * G_G_4_G_N
        return chioverdelta
    

    def get_mu_Sigma_gamma(self):
        """
        Computes deviations in the Poisson and Weyl potentials. Note: current implementation assumes
        G4 = constant, where Sigma = mu.
        """

        if self.output['E'] is None:

            mu, Sigma, gamma = None, None, None

        else:

            if self.output['E'].ndim == 1:

                check = True

                if self.output['beta'] is None or np.isfinite(self.output['beta']).all() == False:
                    check = False

                if check:
                    mu = 1 + self.output['beta']
                    if self.symfunc['alpha_M'].is_zero:
                        Sigma = np.copy(mu)
                        gamma = np.ones(len(mu))
                    else:
                        # using equation 6 and 7 from https://arxiv.org/pdf/2401.06221
                        beta1 = 3*(self.output['rho_phi']/(self.output['H']**2))*(1 + self.output['w_phi']) + (self.output['E_prime']/self.output['E'])*(4 + self.output['alpha_B']) + self.output['alpha_B_prime']
                        beta2 = self.output['alpha_B'] + 2*self.output['alpha_M']
                        factor = 2*beta1 + 2*(1+self.output['alpha_M'])*beta2
                        factor /= 2*beta1 + (2+self.output['alpha_M'])*beta2
                        gamma = 2./factor - 1.
                        Sigma = factor*np.copy(mu)
                else:
                    mu, Sigma, gamma = None, None, None

            elif self.output['E'].ndim == 2:

                mu = np.zeros(np.shape(self.output['E']))
                gamma = np.zeros(np.shape(self.output['E']))
                Sigma = np.zeros(np.shape(self.output['E']))

                for idx in range(0, len(self.output['Omega_c0'])):

                    check = True

                    if self.output['beta'][idx] is None or np.isfinite(self.output['beta'][idx]).all() == False:
                        check = False

                    if check:
                        mu[idx] = 1 + self.output['beta'][idx]
                        if self.symfunc['alpha_M'].is_zero:
                            Sigma[idx] = np.copy(mu[idx])
                            gamma[idx] = np.ones(len(mu[idx]))
                        else:
                            # using equation 6 and 7 from https://arxiv.org/pdf/2401.06221
                            beta1 = 3*(self.output['rho_phi'][idx]/(self.output['H'][idx]**2))*(1 + self.output['w_phi'][idx]) + (self.output['E_prime'][idx]/self.output['E'][idx])*(4 + self.output['alpha_B'][idx]) + self.output['alpha_B_prime'][idx]
                            beta2 = self.output['alpha_B'][idx] + 2*self.output['alpha_M'][idx]
                            factor = 2*beta1 + 2*(1+self.output['alpha_M'][idx])*beta2
                            factor /= 2*beta1 + (2+self.output['alpha_M'][idx])*beta2
                            gamma[idx] = 2./factor - 1.
                            Sigma[idx] = factor*np.copy(mu[idx])
                    else:
                        mu[idx] = np.nan * np.ones(len(self.output['x']))
                        gamma[idx] = np.nan * np.ones(len(self.output['x']))
                        Sigma[idx] = np.nan * np.ones(len(self.output['x']))
            
            else:
                mu, Sigma, gamma = None, None, None
        
        self.output['mu'] = mu
        self.output['gamma'] = gamma
        self.output['Sigma'] = Sigma
    
    # The main solver function

    def _run_solver_start(self, z_max, Npoints, forwards):
        """
        Initialises the solver redshift, scale factor and log(a) x-axes and obtains initial guesses for Hubble and Omegas.

        z_max : float
            Maximum redshift.
        Npoints : int
            Number of points to evaluate numerical functions, from zmax to redshift 0.
        forwards : bool
            Defines whether the solver runs forwards in time (high redshift to low) or backwards.
        """
        # defining redshift range
        z_min = 0.
        x_max = redshift.z2x(z_max)
        x_min = redshift.z2x(z_min)
        x_arr = np.linspace(x_min, x_max, Npoints)
        a_arr = redshift.x2a(x_arr)
        z_arr = redshift.a2z(a_arr)

        if forwards:
            x_arr = x_arr[::-1]
            a_arr = a_arr[::-1]
            z_arr = z_arr[::-1]

        self.output['x'] = x_arr
        self.output['a'] = a_arr
        self.output['z'] = z_arr

        a_start = a_arr[0]
        z_start = z_arr[0]

        self.output['initialiser'] = {}
        self.output['initialiser']['z_start'] = z_start
        self.output['initialiser']['forwards'] = forwards

        # Let's guess the values of the variables by assuming the solution lies close to the reference LCDM values.
        E_ini = self.compute_E_LCDM(a_start)
        E_prime_ini = self.compute_E_prime_LCDM(a_start)

        rho_g_ini = self.get_rho_g(a_start, self.params['H0_ref']*1e-2)
        rho_b_ini = self.get_rho_b(a_start, self.params['Omega_b0_ref'])
        rho_c_ini = self.get_rho_c(a_start, self.params['Omega_c0_ref'])
        rho_l_ini = self.get_rho_l(a_start, self.params['Omega_l0_ref'], w0=self.params['w0'], wa=self.params['wa'])
        rho_nu_ur_ini = self.get_rho_nu_ur(a_start, self.params['H0_ref']*1e-2)
        rho_nu_nr_ini = 0.
        for _mnu in self.params['mnu']:
            rho_nu_nr_ini += self.get_rho_nu_nr(a_start, self.params['H0_ref']*1e-2, _mnu)

        w_nu_nr_ini = 0.
        for _mnu in self.params['mnu']:
            w_nu_nr_ini += self.compute_w_nu_nr(a_start, _mnu)*self.get_rho_nu_nr(a_start, self.params['H0_ref']*1e-2, _mnu)
        if rho_nu_nr_ini != 0:
            w_nu_nr_ini /= rho_nu_nr_ini
        
        w_l_ini   = self.compute_w_l(a_start)
        
        return E_ini, E_prime_ini, rho_g_ini, rho_b_ini, rho_c_ini, rho_l_ini, w_l_ini, rho_nu_ur_ini, rho_nu_nr_ini, w_nu_nr_ini
    

    def _run_solver_initialiser(
            self, 
            variable1, 
            E_ini, E_prime_ini, phi_ini, phi_prime_ini, 
            rho_g_ini, rho_b_ini, rho_c_ini, rho_l_ini, w_l_ini,
            rho_nu_ur_ini, rho_nu_nr_ini, w_nu_nr_ini,
            variable2=None
        ):
        """
        Finds the real roots for initialising the solver for 1 variable, using only the closure, 
        or 2 variables jointly fitting to closure and the E_prime equation.

        Parameters
        ----------
        variable1 : int
            Variable used to set the initial conditions.
        E_ini : float
            Initial normalised Hubble factor.
        E_prime_ini : float
            The derivative of the initial normalised Hubble factor, only used to jointly solve for variable1
            and variable2, if variable2 is None, then E_prime_ini is ignored and solved directly.
        phi_ini : float
            phi initial value, if variable1/variable2 = 0 then this is ignored and solved.
        phi_prime_ini : float
            phi_prime initial value, if variable1/variable2 = 1 then this is ignored and solved.
        rho_g_ini : float
            Initial fractional photon density.
        rho_b_ini : float
            Initial fractional baryon density.
        rho_c_ini : float
            Initial fractional cold dark matter density.
        rho_l_ini : float
            Initial fractional dark energy density.
        w_l_ini : float
            Initial dark energy equation of state.
        rho_nu_ur_ini : float
            Initial density for massless neutrinos.
        rho_nu_nr_ini : float
            Initial density for massive neutrinos.
        w_nu_nr_ini : float
            Initial equation of state for massive neutrinos.
        variable2 : int
            Second variable jointly solved via the closure and E_prime equation to set up the initial conditions.
        """

        # Basic checks

        if variable2 is not None:
            assert variable1 != variable2, "variable1 and 2 must not be the same."
            
        assert variable1 >= 0 and variable1 <= 1, "Variable1 unsupported, must be between 0 and 1 inclusive."
        if variable2 is not None:
            assert variable2 >= 0 and variable2 <= 1, "Variable2 unsupported, must be between 0 and 1 inclusive."

        sub_dict = {}
        
        if variable1 != 0 and variable2 != 0:
            sub_dict[self.sym['phi']] = phi_ini
        else:
            if variable1 == 0:
                variable1_str = 'phi'
                variable1_sym = self.sym['phi']
                variable1_ini = phi_ini
            elif variable2 == 0:
                variable2_str = 'phi'
                variable2_sym = self.sym['phi']
                variable2_ini = phi_ini

        if variable1 != 1 and variable2 != 1:
            sub_dict[self.sym['phi_prime']] = phi_prime_ini
        else:
            if variable1 == 1:
                variable1_str = 'phi_prime'
                variable1_sym = self.sym['phi_prime']
                variable1_ini = phi_prime_ini
            elif variable2 == 1:
                variable2_str = 'phi_prime'
                variable2_sym = self.sym['phi_prime']
                variable2_ini = phi_prime_ini
        
        sub_dict[self.sym['E']] = E_ini
        sub_dict[self.sym['rho_g']] = rho_g_ini
        sub_dict[self.sym['rho_b']] = rho_b_ini
        sub_dict[self.sym['rho_c']] = rho_c_ini
        sub_dict[self.sym['rho_l']] = rho_l_ini
        sub_dict[self.sym['w_l']] = w_l_ini
        sub_dict[self.sym['rho_n_ur']] = rho_nu_ur_ini
        sub_dict[self.sym['rho_n_nr']] = rho_nu_nr_ini
        sub_dict[self.sym['w_n_nr']] = w_nu_nr_ini
        
        self.params['f_H_value'] = 1.
        sub_dict[self.sym['f_H']] = self.params['f_H_value']

        for i in range(len(self.sym['K_G3_G4_syms'])):
            sub_dict[self.sym['K_G3_G4_syms'][i]] = self.params['K_G3_G4_values'][i]

        if variable2 is None:

            closure_func = sym.simplify(self.symfunc['fried_closure'].subs(sub_dict))
            roots1_raw = sym.solve(sym.Eq(closure_func, 0), variable1_sym)

            if self.verbose:
                print(' -- Closure solution for %s_ini:' % variable1_str, roots1_raw)

            roots1 = []
            for root in roots1_raw:
                if root.is_real:
                    roots1.append(root)
            
            variable2_str = None
            roots2_raw = None
            roots2 = None

            if self.verbose:
                print(' -- Real roots for %s_ini:' % variable1_str, roots1)
            
            if len(roots1) > 0:
                ini_success = True
            else:
                ini_success = False

        else:
            
            from scipy.optimize import least_squares

            closure_func = sym.simplify(self.symfunc['fried_closure'].subs(sub_dict))
            E_prime_func = sym.simplify(self.symfunc['E_prime'].subs(sub_dict))

            closure_lambda = sym.lambdify([variable1_sym, variable2_sym], closure_func, 'numpy')
            E_prime_lambda = sym.lambdify([variable1_sym, variable2_sym], E_prime_func, 'numpy')

            def _residuals(variables, E_prime):
                """
                Computes residuals for closure and E_prime equations.

                Parameters
                ----------
                variables : array
                    Variable values
                E_prime : float
                    E_prime to be fitted.
                """
                res0 = closure_lambda(*variables)
                res1 = E_prime_lambda(*variables) - E_prime
                return np.array([res0, res1])
            
            variable0 = [variable1_ini, variable2_ini]   # initial guess
            result = least_squares(_residuals, variable0, args=(E_prime_ini,))

            if result.success:
                ini_success = True
                roots1 = [result.x[0]]
                roots2 = [result.x[1]]
            else:
                ini_success = False
                roots1 = []
                roots2 = []
            
            roots1_raw = roots1
            roots2_raw = roots2

            if self.verbose:
                print(' -- Real roots for %s_ini:' % variable1_str, roots1)
                print(' -- Real roots for %s_ini:' % variable2_str, roots2)
            
        self.output['initialiser']['variable1'] = variable1
        self.output['initialiser']['variable1_string'] = variable1_str
        self.output['initialiser']['all_roots1'] = roots1_raw
        self.output['initialiser']['roots1'] = roots1

        self.output['initialiser']['variable2'] = variable2
        self.output['initialiser']['variable2_string'] = variable2_str
        self.output['initialiser']['all_roots2'] = roots2_raw
        self.output['initialiser']['roots2'] = roots2

        self.output['initialiser']['success'] = ini_success
        self.output['success'] = ini_success
        
        self.output['K_G3_G4_variables'] = [str(sym) for sym in self.sym['K_G3_G4_syms']],
        self.output['K_G3_G4_values'] = self.params['K_G3_G4_values']


    def _run_solver_ODE_NONE(self):
        """
        Returns None for solver outputs.
        """
        self.output['solver_success'] = False
        self.output['success'] = False
        self.output['H0'] = None
        self.output['fH'] = None
        self.output['Omega_g0'] = None
        self.output['Omega_b0'] = None
        self.output['Omega_c0'] = None
        self.output['Omega_l0'] = None
        self.output['Omega_nu_ur0'] = None
        self.output['Omega_nu_nr0'] = None
        self.output['Ehat'] = None
        self.output['Ehat_prime'] = None
        self.output['phihat'] = None
        self.output['phihat_prime'] = None
        self.output['phihat_primeprime'] = None
        self.output['E'] = None
        self.output['E_prime'] = None
        self.output['phi'] = None
        self.output['phi_prime'] = None
        self.output['phi_primeprime'] = None
        self.output['rhohat_phi'] = None
        self.output['rhohat_g'] = None
        self.output['rhohat_b'] = None
        self.output['rhohat_c'] = None
        self.output['rhohat_l'] = None
        self.output['rhohat_nu_ur'] = None
        self.output['rhohat_nu_nr'] = None
        self.output['rho_phi'] = None
        self.output['rho_g'] = None
        self.output['rho_b'] = None
        self.output['rho_c'] = None
        self.output['rho_l'] = None
        self.output['rho_nu_ur'] = None
        self.output['rho_nu_nr'] = None
        self.output['Omega_phi'] = None
        self.output['Omega_g'] = None
        self.output['Omega_b'] = None
        self.output['Omega_c'] = None
        self.output['Omega_l'] = None
        self.output['Omega_nu_ur'] = None
        self.output['Omega_nu_nr'] = None
        self.output['w_nu_nr'] = None
        self.output['w_l'] = None
    

    def _run_solver_ODE_HG(
            self, phi_ini, phi_prime_ini, method='RK45', timeout=5, store_hat=False
        ):
        """
        Returns the Horndeski solver outputs.

        Parameters
        ----------
        phi_ini : float
            phi initial value.
        phi_prime_ini : float
            phi_prime initial value.
        method : str, optional
            solve_ivp method for numerical integration, use 'RK45' for general settings but switch to 'LSODA' if the solver hangs.
        timeout : float, optional
            Time in seconds to force the solver to exit and return nan, this has been added to prvent `solve_ivp` from hanging due 
            to certain variables approaching infinity. You can use different solvers, see method keyword arguement, but this will
            only work if the reason for the failure is due to the equations becoming stiff.
        store_hat : bool, optional
            If true will store raw ODE outputs before normalisation corrections for E renormalisation via f_H.
        """
        if self.output['success'] == False:

            self._run_solver_ODE_NONE()
        
        else:

            x_arr = self.output['x']
            a_arr = self.output['a']
            a_start = self.output['a'][0]
            z_start = self.output['z'][0]

            x_start = x_arr[0]
            x_final = x_arr[-1]

            roots1 = self.output['initialiser']['roots1']
            variable1 = self.output['initialiser']['variable1']

            roots2 = self.output['initialiser']['roots2']
            variable2 = self.output['initialiser']['variable2']
            
            solver_success = [False for r in roots1]
            Ehat_arr = np.zeros((len(roots1), len(x_arr)))
            Ehat_prime_arr = np.zeros((len(roots1), len(x_arr)))
            phihat_arr = np.zeros((len(roots1), len(x_arr)))
            phihat_prime_arr = np.zeros((len(roots1), len(x_arr)))
            phihat_primeprime_arr = np.zeros((len(roots1), len(x_arr)))
            E_arr = np.zeros((len(roots1), len(x_arr)))
            E_prime_arr = np.zeros((len(roots1), len(x_arr)))
            phi_arr = np.zeros((len(roots1), len(x_arr)))
            phi_prime_arr = np.zeros((len(roots1), len(x_arr)))
            phi_primeprime_arr = np.zeros((len(roots1), len(x_arr)))

            Omega_phi_arr = np.zeros((len(roots1), len(x_arr)))
            Omega_g_arr = np.zeros((len(roots1), len(x_arr)))
            Omega_b_arr = np.zeros((len(roots1), len(x_arr)))
            Omega_c_arr = np.zeros((len(roots1), len(x_arr)))
            Omega_l_arr = np.zeros((len(roots1), len(x_arr)))
            Omega_nu_ur_arr = np.zeros((len(roots1), len(x_arr)))
            Omega_nu_nr_arr = np.zeros((len(roots1), len(x_arr)))
            
            rhohat_phi_arr = np.zeros((len(roots1), len(x_arr)))
            rhohat_g_arr = np.zeros((len(roots1), len(x_arr)))
            rhohat_b_arr = np.zeros((len(roots1), len(x_arr)))
            rhohat_c_arr = np.zeros((len(roots1), len(x_arr)))
            rhohat_l_arr = np.zeros((len(roots1), len(x_arr)))
            rhohat_nu_ur_arr = np.zeros((len(roots1), len(x_arr)))
            rhohat_nu_nr_arr = np.zeros((len(roots1), len(x_arr)))

            rho_phi_arr = np.zeros((len(roots1), len(x_arr)))
            rho_g_arr = np.zeros((len(roots1), len(x_arr)))
            rho_b_arr = np.zeros((len(roots1), len(x_arr)))
            rho_c_arr = np.zeros((len(roots1), len(x_arr)))
            rho_l_arr = np.zeros((len(roots1), len(x_arr)))
            rho_nu_ur_arr = np.zeros((len(roots1), len(x_arr)))
            rho_nu_nr_arr = np.zeros((len(roots1), len(x_arr)))
            w_l_arr = np.zeros((len(roots1), len(x_arr)))
            w_nu_nr_arr = np.zeros((len(roots1), len(x_arr)))
            
            H0 = np.zeros(len(roots1)) 
            Omega_phi0 = np.zeros(len(roots1)) 
            Omega_g0 = np.zeros(len(roots1)) 
            Omega_b0 = np.zeros(len(roots1)) 
            Omega_c0 = np.zeros(len(roots1)) 
            Omega_l0 = np.zeros(len(roots1))
            Omega_nu_ur0 = np.zeros(len(roots1)) 
            Omega_nu_nr0 = np.zeros(len(roots1)) 
            fH = np.ones(len(roots1))

            for idx in range(0, len(roots1)):
                
                if variable1 == 0:
                    phi_ini = roots1[idx]
                elif variable1 == 1:
                    phi_prime_ini = roots1[idx]

                if variable2 == 0:
                    phi_ini = roots2[idx]
                elif variable2 == 1:
                    phi_prime_ini = roots2[idx]

                x_ini = x_start

                from scipy.integrate import solve_ivp

                Y_ini = [phi_ini, phi_prime_ini]

                self._initiate_solver_status()
                self._start_timer()
                
                solution = solve_ivp(
                    self._compute_primes, [x_ini, x_final], Y_ini, t_eval=x_arr, method=method, 
                    args=(timeout,), 
                    rtol = 1e-8
                )
                
                solver_success[idx] = self._solver_success
                
                solution = solution["y"].T
                _phi_arr = solution[:,0]
                _phi_prime_arr = solution[:,1]

                _rho_g_arr = np.zeros_like(_phi_arr)
                _rho_b_arr = np.zeros_like(_phi_arr)
                _rho_c_arr = np.zeros_like(_phi_arr)
                _rho_l_arr = np.zeros_like(_phi_arr)
                _rho_nu_ur_arr = np.zeros_like(_phi_arr)
                _rho_nu_nr_arr = np.zeros_like(_phi_arr)
                
                _w_nu_nr_arr = np.zeros_like(_phi_arr)
                _w_l_arr = np.zeros_like(_phi_arr)

                _E_arr = np.zeros_like(_phi_arr)
                
                for i in range(0, len(_phi_arr)):
                    
                    a = np.exp(x_arr[i])

                    _phi = _phi_arr[i]
                    _phi_prime = _phi_prime_arr[i]

                    _rho_g = self.get_rho_g(a, self.params['H0_ref']*1e-2)
                    _rho_b = self.get_rho_b(a, self.params['Omega_b0_ref'])
                    _rho_c = self.get_rho_c(a, self.params['Omega_c0_ref'])
                    _rho_l = self.get_rho_l(a, self.params['Omega_l0_ref'], w0=self.params['w0'], wa=self.params['wa'])
                    _rho_nu_ur = self.get_rho_nu_ur(a, self.params['H0_ref']*1e-2)
                    _rho_nu_nr = 0.
                    for _mnu in self.params['mnu']:
                        _rho_nu_nr += self.get_rho_nu_nr(a, self.params['H0_ref']*1e-2, _mnu)

                    _w_nu_nr = 0.
                    for _mnu in self.params['mnu']:
                        _w_nu_nr += self.compute_w_nu_nr(a, _mnu)*self.get_rho_nu_nr(a, self.params['H0_ref']*1e-2, _mnu)
                    if _rho_nu_nr != 0:
                        _w_nu_nr /= _rho_nu_nr
                    _w_l = self.compute_w_l(a)
                    
                    variables = [
                        _phi, _phi_prime, _rho_g, _rho_b, _rho_c, _rho_l, _rho_nu_ur, _rho_nu_nr, _w_nu_nr,  _w_l, *self.params['K_G3_G4_values'], self.params['fH']
                    ]
                    
                    _E_guess = self.compute_E_LCDM(a)

                    _E_arr[i] = newton(
                        lambda _E: self.lambda_funcs['fried_closure'](_E, *variables), 
                        _E_guess,
                        fprime = lambda _E: self.lambda_funcs['fried_closure_dE'](_E, *variables), 
                        fprime2 = lambda _E: self.lambda_funcs['fried_closure_dE2'](_E, *variables), 
                        tol=self.newton_tol
                    )

                    _rho_g_arr[i] = _rho_g
                    _rho_b_arr[i] = _rho_b
                    _rho_c_arr[i] = _rho_c
                    _rho_l_arr[i] = _rho_l
                    _rho_nu_ur_arr[i] = _rho_nu_ur
                    _rho_nu_nr_arr[i] = _rho_nu_nr

                    _w_nu_nr_arr[i] = _w_nu_nr
                    _w_l_arr[i] = _w_l

                variables = [
                    _E_arr, _phi_arr, _phi_prime_arr, 
                    _rho_g_arr, _rho_b_arr, _rho_c_arr, _rho_l_arr,
                    _rho_nu_ur_arr, _rho_nu_nr_arr, _w_nu_nr_arr,
                    _w_l_arr, *self.params['K_G3_G4_values'], self.params['fH']
                ]

                _rho_phi_arr = self.lambda_funcs['rho_phi'](*variables)
                _E_prime_arr = self.lambda_funcs['E_prime'](*variables)
                _phi_primeprime_arr = self.lambda_funcs['phi_primeprime'](*variables)

                if len(_E_arr) != len(x_arr):
                    self._solver_success = False

                if self._solver_success == False:
                    split = len(_phi_prime_arr)
                    
                if self._solver_success:
                    Ehat_arr[idx] = _E_arr
                    Ehat_prime_arr[idx] = _E_prime_arr
                    phihat_arr[idx] = _phi_arr
                    phihat_prime_arr[idx] = _phi_prime_arr
                    phihat_primeprime_arr[idx] = _phi_primeprime_arr
                    rhohat_phi_arr[idx] = _rho_phi_arr
                    rhohat_g_arr[idx] = _rho_g_arr
                    rhohat_b_arr[idx] = _rho_b_arr
                    rhohat_c_arr[idx] = _rho_c_arr
                    rhohat_l_arr[idx] = _rho_l_arr
                    rhohat_nu_ur_arr[idx] = _rho_nu_ur_arr
                    rhohat_nu_nr_arr[idx] = _rho_nu_nr_arr
                    Omega_phi_arr[idx] = _rho_phi_arr/(_E_arr**2)
                    Omega_g_arr[idx] = _rho_g_arr/(_E_arr**2)
                    Omega_b_arr[idx] = _rho_b_arr/(_E_arr**2)
                    Omega_c_arr[idx] = _rho_c_arr/(_E_arr**2)
                    Omega_l_arr[idx] = _rho_l_arr/(_E_arr**2)
                    Omega_nu_ur_arr[idx] = _rho_nu_ur_arr/(_E_arr**2)
                    Omega_nu_nr_arr[idx] = _rho_nu_nr_arr/(_E_arr**2)
                    w_nu_nr_arr[idx] = _w_nu_nr_arr
                    w_l_arr[idx] = _w_l_arr
                else:
                    Ehat_arr[idx][:split] = _E_arr
                    Ehat_arr[idx][split:] = np.nan
                    Ehat_prime_arr[idx][:split] = _E_prime_arr
                    Ehat_prime_arr[idx][split:] = np.nan
                    phihat_arr[idx][:split] = _phi_arr
                    phihat_arr[idx][split:] = np.nan
                    phihat_prime_arr[idx][:split] = _phi_prime_arr
                    phihat_prime_arr[idx][split:] = np.nan
                    phihat_primeprime_arr[idx][:split] = _phi_primeprime_arr
                    phihat_primeprime_arr[idx][split:] = np.nan
                    rhohat_phi_arr[idx][:split] = _rho_phi_arr
                    rhohat_phi_arr[idx][split:] = np.nan
                    rhohat_g_arr[idx][:split] = _rho_g_arr
                    rhohat_g_arr[idx][split:] = np.nan
                    rhohat_b_arr[idx][:split] = _rho_b_arr
                    rhohat_b_arr[idx][split:] = np.nan
                    rhohat_c_arr[idx][:split] = _rho_c_arr
                    rhohat_c_arr[idx][split:] = np.nan
                    rhohat_l_arr[idx][:split] = _rho_l_arr
                    rhohat_l_arr[idx][split:] = np.nan
                    rhohat_nu_ur_arr[idx][:split] = _rho_nu_ur_arr
                    rhohat_nu_ur_arr[idx][split:] = np.nan
                    rhohat_nu_nr_arr[idx][:split] = _rho_nu_nr_arr
                    rhohat_nu_nr_arr[idx][split:] = np.nan
                    Omega_phi_arr[idx][:split] = _rho_phi_arr/(_E_arr**2)
                    Omega_phi_arr[idx][split:] = np.nan
                    Omega_g_arr[idx][:split] = _rho_g_arr/(_E_arr**2)
                    Omega_g_arr[idx][split:] = np.nan
                    Omega_b_arr[idx][:split] = _rho_b_arr/(_E_arr**2)
                    Omega_b_arr[idx][split:] = np.nan
                    Omega_c_arr[idx][:split] = _rho_c_arr/(_E_arr**2)
                    Omega_c_arr[idx][split:] = np.nan
                    Omega_l_arr[idx][:split] = _rho_l_arr/(_E_arr**2)
                    Omega_l_arr[idx][split:] = np.nan
                    Omega_nu_ur_arr[idx][:split] = _rho_nu_ur_arr/(_E_arr**2)
                    Omega_nu_ur_arr[idx][split:] = np.nan
                    Omega_nu_nr_arr[idx][:split] = _rho_nu_nr_arr/(_E_arr**2)
                    Omega_nu_nr_arr[idx][split:] = np.nan
                    w_nu_nr_arr[idx][:split] = _w_nu_nr_arr
                    w_nu_nr_arr[idx][split:] = np.nan
                    w_l_arr[idx][:split] = _w_l_arr
                    w_l_arr[idx][split:] = np.nan
                
                if z_start == 0.:
                    self.params['H0'] = self.params['H0_ref']*Ehat_arr[idx][0]
                else:
                    if self._solver_success == True:
                        self.params['H0'] = self.params['H0_ref']*Ehat_arr[idx][-1]
                    else:
                        self.params['H0'] = np.nan
                
                if np.isfinite(self.params['H0']) and self.params['H0'] != self.params['H0_ref'] and self._solver_success:
                    self.params['f_H_value'] = self.params['H0']/self.params['H0_ref']
                else:
                    self.params['f_H_value'] = 1.

                E_arr[idx] = Ehat_arr[idx]/self.params['f_H_value']
                E_prime_arr[idx] = Ehat_prime_arr[idx]/self.params['f_H_value']
                phi_arr[idx] = phihat_arr[idx]*self.params['f_H_value']
                phi_prime_arr[idx] = phihat_prime_arr[idx]*self.params['f_H_value']
                phi_primeprime_arr[idx] = phihat_primeprime_arr[idx]*self.params['f_H_value']
                rho_phi_arr[idx] = rhohat_phi_arr[idx]/(self.params['f_H_value']**2)
                rho_g_arr[idx] = rhohat_g_arr[idx]/(self.params['f_H_value']**2)
                rho_c_arr[idx] = rhohat_b_arr[idx]/(self.params['f_H_value']**2)
                rho_b_arr[idx] = rhohat_c_arr[idx]/(self.params['f_H_value']**2)
                rho_l_arr[idx] = rhohat_l_arr[idx]/(self.params['f_H_value']**2)
                rho_nu_ur_arr[idx] = rhohat_nu_ur_arr[idx]/(self.params['f_H_value']**2)
                rho_nu_nr_arr[idx] = rhohat_nu_nr_arr[idx]/(self.params['f_H_value']**2)

                fH[idx] = self.params['f_H_value']
                H0[idx] = self.params['H0']

                if z_start == 0.:
                    self.params['Omega_phi0'] = Omega_phi_arr[idx][0]
                    self.params['Omega_g0'] = Omega_g_arr[idx][0]
                    self.params['Omega_b0'] = Omega_b_arr[idx][0]
                    self.params['Omega_c0'] = Omega_c_arr[idx][0]
                    self.params['Omega_l0'] = Omega_l_arr[idx][0]
                    self.params['Omega_nu_ur0'] = Omega_nu_ur_arr[idx][0]
                    self.params['Omega_nu_nr0'] = Omega_nu_nr_arr[idx][0]
                else:
                    if self._solver_success == True:
                        self.params['Omega_phi0'] = Omega_phi_arr[idx][-1]
                        self.params['Omega_g0'] = Omega_g_arr[idx][-1]
                        self.params['Omega_b0'] = Omega_b_arr[idx][-1]
                        self.params['Omega_c0'] = Omega_c_arr[idx][-1]
                        self.params['Omega_l0'] = Omega_l_arr[idx][-1]
                        self.params['Omega_nu_ur0'] = Omega_nu_ur_arr[idx][-1]
                        self.params['Omega_nu_nr0'] = Omega_nu_nr_arr[idx][-1]
                    else:
                        self.params['Omega_phi0'] = np.nan
                        self.params['Omega_g0'] = np.nan
                        self.params['Omega_b0'] = np.nan
                        self.params['Omega_c0'] = np.nan
                        self.params['Omega_l0'] = np.nan
                        self.params['Omega_nu_ur0'] = np.nan
                        self.params['Omega_nu_nr0'] = np.nan
                
                Omega_phi0[idx] = self.params['Omega_phi0']
                Omega_g0[idx] = self.params['Omega_g0']
                Omega_b0[idx] = self.params['Omega_b0']
                Omega_c0[idx] = self.params['Omega_c0']
                Omega_l0[idx] = self.params['Omega_l0']
                Omega_nu_ur0[idx] = self.params['Omega_nu_ur0']
                Omega_nu_nr0[idx] = self.params['Omega_nu_nr0']
            
            self.output['solver_success'] = solver_success
            self.output['success'] = solver_success
            self.output['H0'] = H0
            self.output['fH'] = fH
            self.output['Omega_phi0'] = Omega_phi0
            self.output['Omega_g0'] = Omega_g0
            self.output['Omega_b0'] = Omega_b0
            self.output['Omega_c0'] = Omega_c0
            self.output['Omega_l0'] = Omega_l0
            self.output['Omega_nu_ur0'] = Omega_nu_ur0
            self.output['Omega_nu_nr0'] = Omega_nu_nr0
            if store_hat:
                self.output['Ehat'] = Ehat_arr
                self.output['Ehat_prime'] = Ehat_prime_arr
                self.output['phihat'] = phihat_arr
                self.output['phihat_prime'] = phihat_prime_arr
                self.output['phihat_primeprime'] = phihat_primeprime_arr
            else:
                self.output['Ehat'] = None
                self.output['Ehat_prime'] = None
                self.output['phihat'] = None
                self.output['phihat_prime'] = None
                self.output['phihat_primeprime'] = None
            self.output['E'] = E_arr
            self.output['E_prime'] = E_prime_arr
            self.output['phi'] = phi_arr
            self.output['phi_prime'] = phi_prime_arr
            self.output['phi_primeprime'] = phi_primeprime_arr
            if store_hat:
                self.output['rhohat_phi'] = rhohat_phi_arr
                self.output['rhohat_g'] = rhohat_g_arr
                self.output['rhohat_b'] = rhohat_b_arr
                self.output['rhohat_c'] = rhohat_c_arr
                self.output['rhohat_l'] = rhohat_l_arr
                self.output['rhohat_nu_ur'] = rhohat_nu_ur_arr
                self.output['rhohat_nu_nr'] = rhohat_nu_nr_arr
            else:
                self.output['rhohat_phi'] = None
                self.output['rhohat_g'] = None
                self.output['rhohat_b'] = None
                self.output['rhohat_c'] = None
                self.output['rhohat_l'] = None
                self.output['rhohat_nu_ur'] = None
                self.output['rhohat_nu_nr'] = None
            self.output['rho_phi'] = rho_phi_arr
            self.output['rho_g'] = rho_g_arr
            self.output['rho_b'] = rho_b_arr
            self.output['rho_c'] = rho_c_arr
            self.output['rho_l'] = rho_l_arr
            self.output['rho_nu_ur'] = rho_nu_ur_arr
            self.output['rho_nu_nr'] = rho_nu_nr_arr
            self.output['Omega_phi'] = Omega_phi_arr
            self.output['Omega_g'] = Omega_g_arr
            self.output['Omega_b'] = Omega_b_arr
            self.output['Omega_c'] = Omega_c_arr
            self.output['Omega_l'] = Omega_l_arr
            self.output['Omega_nu_ur'] = Omega_nu_ur_arr
            self.output['Omega_nu_nr'] = Omega_nu_nr_arr
            self.output['w_nu_nr'] = w_nu_nr_arr
            self.output['w_l'] = w_l_arr
    

    def _run_solver_derived_NONE(self):
        """
        Returns None for derived outputs.
        """
        self.output['H'] = None
        self.output['Dc'] = None
        self.output['G_G_4/G_N'] = None
        self.output['A'] = None
        self.output['Omega_DE'] = None
        self.output['w_DE'] = None
        self.output['Omega_phi_via_closure'] = None
        self.output['calB'] = None
        self.output['calC'] = None
        self.output['beta'] = None
        self.output['chi/delta'] = None
        self.output['M_star_sq'] = None
        self.output['alpha_M'] = None
        self.output['alpha_B'] = None
        self.output['alpha_B_prime'] = None
        self.output['alpha_K'] = None
        self.output['tilde_calE'] = None
        self.output['tilde_calP'] = None
        self.output['w_phi'] = None
        self.output['D'] = None
        self.output['Q_s'] = None
        self.output['c_s_sq_D'] = None
        self.output['c_s_sq'] = None
        self.output['stable'] = None
    

    def _run_solver_derived_HG(self):
        """
        Returns derived outputs for Horndeski gravity.
        """

        if self.output['success'] == False:

            self._run_solver_derived_NONE()
        
        else:
            
            x_arr = self.output['x']
            a_arr = self.output['a']
            E_arr = self.output['E']
            E_prime_arr = self.output['E']
            phi_arr = self.output['phi']
            phi_prime_arr = self.output['phi_prime']
            rho_g_arr = self.output['rho_g']
            rho_b_arr = self.output['rho_b']
            rho_c_arr = self.output['rho_c']
            rho_l_arr = self.output['rho_l']
            rho_nu_ur_arr = self.output['rho_nu_ur']
            rho_nu_nr_arr = self.output['rho_nu_nr']
            Omega_phi_arr = self.output['Omega_phi']
            Omega_g_arr = self.output['Omega_g']
            Omega_b_arr = self.output['Omega_b']
            Omega_c_arr = self.output['Omega_c']
            Omega_l_arr = self.output['Omega_l']
            Omega_nu_ur_arr = self.output['Omega_nu_ur']
            Omega_nu_nr_arr = self.output['Omega_nu_nr']
            w_nu_nr_arr = self.output['w_nu_nr']
            w_l_arr = self.output['w_l']

            roots1 = self.output['initialiser']['roots1']
            
            G_G_4_G_N = np.zeros((len(roots1), len(x_arr)))
            A_arr = np.zeros((len(roots1), len(x_arr)))
            Omega_DE_arr = np.zeros((len(roots1), len(x_arr)))
            w_DE_arr = np.zeros((len(roots1), len(x_arr)))
            Omega_phi_via_closure_arr = np.zeros((len(roots1), len(x_arr)))

            calB_arr = np.zeros((len(roots1), len(x_arr)))
            calC_arr = np.zeros((len(roots1), len(x_arr)))
            beta_arr = np.zeros((len(roots1), len(x_arr)))
            chioverdelta_arr = np.zeros((len(roots1), len(x_arr)))
            
            M_star_sq_arr = np.zeros((len(roots1), len(x_arr)))
            alpha_M_arr = np.zeros((len(roots1), len(x_arr)))
            alpha_B_arr = np.zeros((len(roots1), len(x_arr)))
            alpha_B_prime_arr = np.zeros((len(roots1), len(x_arr)))
            alpha_K_arr = np.zeros((len(roots1), len(x_arr)))

            tilde_calE_arr = np.zeros((len(roots1), len(x_arr)))
            tilde_calP_arr = np.zeros((len(roots1), len(x_arr)))
            w_phi_arr = np.zeros((len(roots1), len(x_arr)))

            D_arr = np.zeros((len(roots1), len(x_arr)))
            Q_s_arr = np.zeros((len(roots1), len(x_arr)))
            c_s_sq_D_arr = np.zeros((len(roots1), len(x_arr)))
            c_s_sq_arr = np.zeros((len(roots1), len(x_arr)))

            stable = [True for r in roots1]

            for (idx, _) in enumerate(roots1):

                self.params['H0'] = self.output['H0'][idx]
                self.params['f_H_value'] = self.output['fH'][idx]

                variables = [
                    E_arr[idx], phi_arr[idx], phi_prime_arr[idx], 
                    rho_g_arr[idx], rho_b_arr[idx],
                    rho_c_arr[idx], rho_l_arr[idx], 
                    rho_nu_ur_arr[idx],
                    rho_nu_nr_arr[idx], w_nu_nr_arr[idx], 
                    w_l_arr[idx], *self.params['K_G3_G4_values'], self.params['f_H_value']
                ]

                G_G_4_G_N[idx] = self.lambda_funcs['G_G_4/G_N'](*variables)
                A_arr[idx] = self.lambda_funcs['A'](*variables)
                
                Omega_DE_arr[idx] = 1. - Omega_g_arr[idx] - Omega_b_arr[idx] - Omega_c_arr[idx] - Omega_nu_ur_arr[idx] - Omega_nu_nr_arr[idx]
                Omega_phi_via_closure_arr[idx] = 1 - Omega_g_arr[idx] - Omega_b_arr[idx] - Omega_c_arr[idx] - Omega_l_arr[idx] - Omega_nu_ur_arr[idx] - Omega_nu_nr_arr[idx]

                calB_arr[idx] = self.lambda_funcs['calB'](*variables)
                calC_arr[idx] = self.lambda_funcs['calC'](*variables)
                beta_arr[idx] = self.lambda_funcs['beta'](*variables)
                chioverdelta_arr[idx] = self.compute_chi_over_delta(a_arr, E_arr[idx], calB_arr[idx], calC_arr[idx], G_G_4_G_N[idx])

                M_star_sq_arr[idx] = self.lambda_funcs['M_star_sq'](*variables)
                alpha_M_arr[idx] = self.lambda_funcs['alpha_M'](*variables)
                alpha_B_arr[idx] = self.lambda_funcs['alpha_B'](*variables)
                alpha_B_prime_arr[idx] = self.lambda_funcs['alpha_B_prime'](*variables)
                alpha_K_arr[idx] = self.lambda_funcs['alpha_K'](*variables)
                
                tilde_calE_arr[idx] = self.lambda_funcs['tilde_calE'](*variables)
                tilde_calP_arr[idx] = self.lambda_funcs['tilde_calP'](*variables)

                # given by eq. 3.10 of Bellini and Sawicki 2014
                w_phi_arr[idx] = tilde_calP_arr[idx]/tilde_calE_arr[idx]
                
                w_DE_arr[idx] = self.compute_w_l(a_arr)*Omega_l_arr[idx] + w_phi_arr[idx]*Omega_phi_arr[idx]
                w_DE_arr[idx] /= Omega_l_arr[idx] + Omega_phi_arr[idx]

                D_arr[idx] = self.lambda_funcs['D'](*variables)
                Q_s_arr[idx] = self.lambda_funcs['Q_s'](*variables)
                c_s_sq_D_arr[idx] = self.lambda_funcs['c_s_sq_D'](*variables)
                c_s_sq_arr[idx] = self.lambda_funcs['c_s_sq'](*variables)

                if Q_s_arr.all() > 0 and c_s_sq_arr.all() > 0:
                    stable[idx] = True
                else:
                    stable[idx] = False

            H_arr = np.zeros_like(E_arr)
            Dc_arr = np.zeros_like(E_arr)
            
            for (idx, _) in enumerate(roots1):
                
                self.params['H0'] = self.output['H0'][idx]
                H_arr[idx] = self.params['H0']*E_arr[idx]

                from scipy.integrate import cumulative_trapezoid
                
                f = 1./E_arr[idx]
                if self.output['initialiser']['forwards']:
                    Dc_arr[idx] = self.const['Dh']*cumulative_trapezoid(f[::-1], x=self.output['z'][::-1], initial=0.)[::-1]
                else:
                    Dc_arr[idx] = self.const['Dh']*cumulative_trapezoid(f, x=self.output['z'], initial=0.)

            self.output['H'] = H_arr
            self.output['Dc'] = Dc_arr
            self.output['G_G_4/G_N'] = G_G_4_G_N
            self.output['A'] = A_arr
            self.output['Omega_DE'] = Omega_DE_arr
            self.output['w_DE'] = w_DE_arr
            self.output['Omega_phi_via_closure'] = Omega_phi_via_closure_arr
            self.output['calB'] = calB_arr
            self.output['calC'] = calC_arr
            self.output['beta'] = beta_arr
            self.output['chi/delta'] = chioverdelta_arr
            self.output['M_star_sq'] = M_star_sq_arr
            self.output['alpha_M'] = alpha_M_arr
            self.output['alpha_B'] = alpha_B_arr
            self.output['alpha_B_prime'] = alpha_B_prime_arr
            self.output['alpha_K'] = alpha_K_arr
            self.output['tilde_calE'] = tilde_calE_arr
            self.output['tilde_calP'] = tilde_calP_arr
            self.output['w_phi'] = w_phi_arr
            self.output['D'] = D_arr
            self.output['Q_s'] = Q_s_arr
            self.output['c_s_sq_D'] = c_s_sq_D_arr
            self.output['c_s_sq'] = c_s_sq_arr
            self.output['stable'] = stable


    def run_solver(
            self, z_max=2000., Npoints=1000, forwards=True, GR=False, variable1=2, variable2=None, 
            phi_ini=1e-6, phi_prime_ini=1e-6, method='RK45', timeout=5, newton_tol=1e-5,
            derived=True, LCDM_ini=True, values_ini=None, store_hat=False, HS_correction=True
        ):
        """
        Runs the numerical solver for a user defined Horndeski model.

        Parameters
        ----------
        z_max : float, optional
            Maximum redshift.
        Npoints : int, optional
            Number of points to evaluate numerical functions, from zmax to redshift 0.
        forwards : bool, optional
            Defines whether the solver runs forwards in time (high redshift to low) or backwards.
        GR : bool, optional
            Force to run with general relativity equations.
        variable1 : str, optional
            Variable used to set the initial conditions, default=1 => phi_prime.
        variable2 : str, optional
            Second variable used to set the initial conditions, fitted jointly with the first default=None.
        phi_ini : float, optional
            phi initial value, if the closure_variable = 1 then this is used as an initial guess.
        phi_prime_ini : float, optional
            phi_prime initial value, if the closure_variable = 2 then this is used as an initial guess.
        method : str, optional
            solve_ivp method for numerical integration, use 'RK45' for general settings but switch to 'LSODA' if the solver hangs.
        timeout : float, optional
            Time in seconds to force the solver to exit and return nan, this has been added to prvent `solve_ivp` from hanging due 
            to certain variables approaching infinity. You can use different solvers, see method keyword arguement, but this will
            only work if the reason for the failure is due to the equations becoming stiff.
        derived : bool, optional
            Instructs the solver whether derived quantities should be computed.
        LCDM_ini : bool, optional
            Instructs the solver to use LCDM initial conditions.
        values_ini : bool, optional
            Directly supply initial values for initial conditions for [E_ini, E_prime_ini, Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini, 
            w_l_ini, Omega_nu_ur, Omega_nu_nr, w_nu_nr_ini].
        store_hat : bool, optional
            If true will store raw ODE outputs before normalisation corrections for E renormalisation via f_H.
        HS_correction : bool, optional
            Applies a bias correction to the Hu & Sugiyama prediction for z_star which is only valid
            for models close to Planck LCDM values during the early universe.
        
        Returns
        -------
        outputs : dict
            Dictionary containing numerical solver solutions.
        """

        if self.verbose:
            print('Hi-COLA: Running numerical ODE solver')
            print(' - Initialising solver...')

        self.output = {}

        (E_ini, E_prime_ini, Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini,
         w_l_ini, Omega_nu_ur_ini, Omega_nu_nr_ini, w_nu_nr_ini) = \
            self._run_solver_start(z_max, Npoints, forwards)
        
        if LCDM_ini == False:
            assert values_ini is not None, (
                "values_ini must be a list with the initial values for "
                "[E_ini, E_prime_ini, Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini, "
                "w_l_ini, Omega_nu_ur, Omega_nu_nr, w_nu_nr_ini]"
            )
            (E_ini, E_prime_ini,
             Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini, w_l_ini,
             Omega_nu_ur_ini, 
             Omega_nu_nr_ini, w_nu_nr_ini) = values_ini
        
        if GR:
            if self.verbose:
                print(' - Running in GR mode!')

            if self.verbose:
                print(' - Computing expansion history.')
            
            self._run_solver_ODE_GR()

            if derived:
                if self.verbose:
                    print(' - Computing main derived quantities.')
                
                self._run_solver_derived_GR()

        else:
            if self.verbose:
                print(' - Running in Horndeski mode!')

            self._run_solver_initialiser(
                variable1,
                E_ini, E_prime_ini, phi_ini, phi_prime_ini, 
                Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini, w_l_ini,
                Omega_nu_ur_ini, Omega_nu_nr_ini, w_nu_nr_ini,
                variable2,
            )

            if self.output['success']:
                
                if self.verbose:
                    print(' - Numerically solving ODEs for expansion history and scalar field evolution.')

                # set Newton-Raphson tolerance for computing E from closure
                self.newton_tol = newton_tol

                self._run_solver_ODE_HG(
                    phi_ini, phi_prime_ini,
                    method=method, timeout=timeout, 
                    store_hat=store_hat
                )

                if derived:
                    if self.verbose:
                        print(' - Computing main derived quantities.')

                    self._run_solver_derived_HG()

            else:
                
                if self.verbose:
                    print(' - No roots found. Solver stopped')
                
                self._run_solver_ODE_NONE()

                if derived:

                    self._run_solver_derived_NONE()
        
        if forwards is False:

            if self.verbose:
                print(' - Reversing ordering for backwards solve.')

            self._reverse_outputs(derived)
        
        if derived:

            if self.verbose:
                print(' - Computing additional derived quantities numerically.')
        
            # compute the effective equation of state
            self.comp_w_eff()

            # compute the phi sound horizon
            self.get_phi_sound_horizon()

            # Compute growth functions
            self.get_linear_growth()
            self.get_linear_growth_2()

            # compute mu, Sigma and gamma variables
            self.get_mu_Sigma_gamma()
            self.get_Sigma_derivatives()

            # compute BA0 related quantities
            self.get_r_drag()
            self.get_a_eq()
            self.get_z_star(apply_correction=HS_correction)
            self.get_r_star()
            self.get_theta_star()


        if self.verbose:
                print(' - Done!')
        
        return self.output


    # Save backend files
    
    def save_backend(self, fname_prefix):
        """
        Save outputs for background Hi-COLA run.

        Parameters
        ----------
        fname_prefix : str
            Filename for Hi-COLA backend files. This will create a force and expansion file in ascii format.
        """
        
        if fname_prefix[-1] != '_':
            fname_prefix += '_'
        
        if utils.isscalar(self.output['Omega_c0']):

            fname_expansion = fname_prefix + 'expansion.txt'
            data = np.column_stack([self.output['a'], self.output['E'], self.output['E_prime']/self.output['E']])
            np.savetxt(fname_expansion, data, fmt=['%.4e', '%.4e', '%.4e'])

            fname_force = fname_prefix + 'force.txt'
            data = np.column_stack([self.output['a'], self.output['chi/delta'], self.output['beta']])
            np.savetxt(fname_force, data, fmt=['%.4e', '%.4e', '%.4e'])
        
        else:
            
            for idx in range(0, len(self.output['Omega_c0'])):
                
                if len(self.output['Omega_c0']) != 1:
                    fname_expansion = fname_prefix + 'solution_%i_expansion.txt' % idx
                else:
                    fname_expansion = fname_prefix + 'expansion.txt'
                data = np.column_stack([self.output['a'], self.output['E'][idx], self.output['E_prime/E'][idx]/self.output['E'][idx]])
                np.savetxt(fname_expansion, data, fmt=['%.4e', '%.4e', '%.4e'])

                if len(self.output['Omega_c0']) != 1:
                    fname_force = fname_prefix + 'solution_%i_force.txt' % idx
                else:
                    fname_force = fname_prefix + 'force.txt'
                data = np.column_stack([self.output['a'], self.output['chi/delta'][idx], self.output['beta'][idx]])
                np.savetxt(fname_force, data, fmt=['%.4e', '%.4e', '%.4e'])


    # Save complete outputs

    def save_outputs(self, fname_prefix, save_backend=True):
        """
        Save output dictionary into numpy format.

        Parameters
        ----------
        fname_prefix : str
            Filename prefix for output file, a 'all.npz' will be added to the name of the file.
        """
        if fname_prefix[-1] != '_':
            fname_prefix += '_'

        fname_all = fname_prefix + 'all.npz'

        np.savez(fname_all, **self.output)
        if save_backend:
            self.save_backend(fname_prefix)

    
    # Reinitialise the class

    def clean(self):
        """
        Reinitialise the class. 
        """
        self.__init__()
