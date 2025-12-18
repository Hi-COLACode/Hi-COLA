import numpy as np
import time
import sympy as sym

from . import lcdm, redshift


class HorndeskiModel:

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
        # Enable print statements
        self.verbose = True
        # Construct dictionary containing sympy symbolic variables.
        self.sym = {}
        self.sym['a'] = sym.symbols('a')
        self.sym['f_H'] = sym.symbols('f_H')
        self.sym['E'] = sym.symbols('E')
        self.sym['E_prime'] = sym.symbols("E^{'}")
        self.sym['phi'] = sym.symbols('phi')
        self.sym['phi_prime'] = sym.symbols("phi^{'}")
        self.sym['phi_primeprime'] = sym.symbols("phi^{''}")
        self.sym['X'] = sym.symbols('X')
        self.sym['M_p'] = sym.symbols('M_{p}')
        self.sym['M_s'] = sym.symbols('M_{s}')
        self.sym['M_g'] = sym.symbols('M_{g}')
        self.sym['M_G3'] = sym.symbols('M_{G3}')
        self.sym['M_G4'] = sym.symbols('M_{G4}')
        self.sym['M_K'] = sym.symbols('M_{K}')
        self.sym['M_pG4'] = sym.symbols('M_{pG4}')
        self.sym['M_KG4'] = sym.symbols('M_{KG4}')
        self.sym['M_G3s'] = sym.symbols('M_{G3s}')
        self.sym['M_sG4'] = sym.symbols('M_{sG4}')
        self.sym['M_G3G4'] = sym.symbols('M_{G3G4}')
        self.sym['M_Ks'] = sym.symbols('M_{Ks}')
        self.sym['M_gp'] = sym.symbols('M_{gp}')
        self.sym['M_sp'] = sym.symbols('M_{sp}')
        self.sym['M_Kp2'] = sym.symbols('M_{Kp2}')
        self.sym['M_G3p'] = sym.symbols('M_{G3p}')
        self.sym['M_G4p2'] = sym.symbols('M_{G4p2}')
        self.sym['Omega_g'] = sym.symbols('Omega_g') # Photons
        self.sym['Omega_n1'] = sym.symbols('Omega_n1') # Neutrino 1
        self.sym['Omega_n2'] = sym.symbols('Omega_n2') # Neutrino 2
        self.sym['Omega_n3'] = sym.symbols('Omega_n3') # Neutrino 3
        self.sym['w_n1'] = sym.symbols('w_n1') # Neutrino 1 EoS
        self.sym['w_n2'] = sym.symbols('w_n2') # Neutrino 2 EoS
        self.sym['w_n3'] = sym.symbols('w_n3') # Neutrino 3 EoS
        self.sym['Omega_b'] = sym.symbols('Omega_b') # Baryon
        self.sym['Omega_c'] = sym.symbols('Omega_c') # Cold dark matter
        self.sym['Omega_l'] = sym.symbols('Omega_l') # Lambda or Dynamical Dark Energy
        self.sym['w_l'] = sym.symbols('w_l') # Lambda or Dynamical Dark Energy EoS
        self.sym['H0'] = sym.symbols('H_0')
        self.sym['f_phi'] = sym.symbols('f_phi')
        self.sym['Theta'] = sym.symbols('Theta')
        # Not sure these are variables or diagnostic tools...
        self.sym['threshold'] = sym.symbols('threshold')
        self.sym['threshold_sign'] = sym.symbols('threshold_sign')
        self.symfunc = {}
        # Lambda functions
        self.lambda_funcs = {}
        # Parameter values
        self.params = {}
        self.params['mass_ratios'] = None
        self.params['mass_ratios_new'] = None
        self.outputs = {}
        self.timer = {'t0': None}
        self.const = {}
        self.const['c[m/s]'] = 299_792_458 # Units m/s
        self.const['c[km/s]'] = 299_792.458 # Units km/s
        self.const['Dh'] = 1e-2 * self.const['c[km/s]'] # Units Mpc/h
        self.const['G'] = 6.6743015e-11 
        self.const['kB'] = 1.380649e-23
        self.const['eV'] = 1.602176634e-19
        self.const['hbar'] = 1.054571817e-34
        self.const['Mpc'] = 3.0857e22
        self.const['kB[eV]'] = self.const['kB']/self.const['eV']
        self.const['m[eV]'] = self.const['eV']/(self.const['c[m/s]']*self.const['hbar']) # meters in eV
        self.const['s[eV]'] = self.const['eV']/self.const['hbar'] # seconds in eV
        self.const['kg[eV]'] = self.const['c[m/s]']**2 / self.const['eV'] # kilograms in eV
        self.const['G[eV]'] = self.const['G'] * (self.const['m[eV]']**3)/(self.const['kg[eV]']*self.const['s[eV]']**2)
        self.const['H0unit'] = 100.*1e3/self.const['Mpc'] # units h.s^-1
        self.const['H0unit[eV]'] = self.const['H0unit']/self.const['s[eV]']
        from scipy.special import zeta
        self.const['riemann_zeta3'] = zeta(3)
        self.const['riemann_zeta5'] = zeta(5)
        # Neutrino table
        self.neutrino_table = None


    # Timer related functions

    def _start_timer(self):
        self.timer['t0'] = time.time()
    

    def _check_timer(self):
        timenow = time.time()
        return timenow - self.timer['t0']
    

    # Clean parameter values

    def clean_params(self):
        """
        Reinitialised the class parameters.
        """
        self.params = {}
        self.outputs = {}
    
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
        self.symfunc['K'] = sym.sympify(K_func)/(self.sym['f_H']**2)
    

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
        self.symfunc['G3'] = sym.sympify(G3_func)/self.sym['f_H']
    
    
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
        # Note: replaces K_func
        self.symfunc['Kx'] = sym.diff(self.symfunc['K'],self.sym['X'])
        self.symfunc['Kphi'] = sym.diff(self.symfunc['K'],self.sym['phi'])
        self.symfunc['Kxx'] = sym.diff(self.symfunc['Kx'],self.sym['X'])
        self.symfunc['Kxphi'] = sym.diff(self.symfunc['Kx'],self.sym['phi'])


    def get_G3_derivatives(self):
        """
        Computes the symbolic differentials of the G3 function with respect to X and phi.
        """
        # Note: replaces G3_func
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
        # Note: replaces G4_func
        self.symfunc['G4x'] = sym.diff(self.symfunc['G4'],self.sym['X'])
        self.symfunc['G4xx'] = sym.diff(self.symfunc['G4x'],self.sym['X'])
        self.symfunc['G4xphi'] = sym.diff(self.symfunc['G4x'],self.sym['phi'])
        self.symfunc['G4phi'] = sym.diff(self.symfunc['G4'],self.sym['phi'])
        self.symfunc['G4phiphi'] = sym.diff(self.symfunc['G4phi'],self.sym['phi'])
        self.symfunc['G4phix'] = sym.diff(self.symfunc['G4phi'],self.sym['X'])

    # Symbolic function definitions

    def get_Omega_phi(self):
        """
        Returns the scalar field fractional density function, following equation 2.4 in https://arxiv.org/pdf/2209.01666.pdf.
        """
        # Note: replaces the omega_phi function
        term1 = (self.sym['Omega_g'] + self.sym['Omega_n1'] + self.sym['Omega_n2'] + self.sym['Omega_n3'] + self.sym['Omega_b'] + self.sym['Omega_c'] + self.sym['Omega_l'])*((self.sym['M_pG4']**2.)/(2.*self.symfunc['G4']) - 1.)
        term2 = (self.sym['M_KG4']**2.)*self.sym['X']*self.symfunc['Kx']/(self.sym['E']**2.) - (self.sym['M_KG4']**2.)*self.symfunc['K']/(2.*(self.sym['E']**2.))
        term2 += 3*self.sym['M_G3G4']*self.sym['M_sG4']*self.sym['X']*self.sym['phi_prime']*self.symfunc['G3x']
        term2 += -self.sym['M_G3G4']*self.sym['M_sG4']*self.sym['X']*self.symfunc['G3phi']/(self.sym['E']**2.) - 3*self.sym['phi_prime']*self.symfunc['G4phi']
        self.symfunc['Omega_phi'] = term1 + (1/(3.*self.symfunc['G4']))*term2

    
    def get_Omega_phi_new(self):
        """
        Returns the scalar field fractional density function, following equation 2.4 in https://arxiv.org/pdf/2209.01666.pdf.
        """
        self.symfunc['Omega_phi_new'] = (1/(2*self.sym['M_G4p2']*self.symfunc['G4']) - 1)
        self.symfunc['Omega_phi_new'] *= (self.sym['Omega_g'] + self.sym['Omega_n1'] + self.sym['Omega_n2'] + self.sym['Omega_n3'] + self.sym['Omega_b'] + self.sym['Omega_c'] + self.sym['Omega_l'])
        term1 = self.sym['M_Kp2']*self.sym['X']*self.symfunc['Kx']/(self.sym['E']**2)
        term1 -= self.sym['M_Kp2']*self.symfunc['K']/(2*self.sym['E']**2)
        term1 += 3*self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.sym['phi_prime']*self.symfunc['G3x']
        term1 -= self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.symfunc['G3phi']/(self.sym['E']**2)
        term1 -= 3*self.sym['M_G4p2']*self.sym['phi_prime']*self.symfunc['G4phi']
        self.symfunc['Omega_phi_new'] += (1/(3*self.sym['M_G4p2']*self.symfunc['G4']))*term1

    def get_fried_closure_new(self):
        """
        Returns the RHS of the Friedmann closure relation -> equation 2.3 in https://arxiv.org/abs/2209.01666, 
        should be equal to 1. 
        """
        # Note: replaces the fried_closure function
        self.get_Omega_phi_new()
        self.symfunc['fried_closure_new'] = self.symfunc['Omega_phi_new'] + self.sym['Omega_b'] + self.sym['Omega_c'] + self.sym['Omega_g'] + self.sym['Omega_n1'] + self.sym['Omega_n2'] + self.sym['Omega_n3'] + self.sym['Omega_l'] - 1.
        # equation A.3 in https://arxiv.org/abs/2209.01666, 
        # TODO: duplication of a step done later in construct model, remove these lines.
        # Xreal = (1./2.)*(self.sym['E']**2.)*(self.sym['phi_prime']**2.)
        # self.symfunc['fried_closure'] = self.symfunc['fried_closure'].subs(self.sym['X'], Xreal)


    def get_G_G_4(self):
        """
        Returns the G_G_4/G_N function.
        """
        self.symfunc['G_G_4/G_N'] = 1/(2*self.symfunc['G4'])
    
    def get_G_G_4_new(self):
        """
        Returns the G_G_4/G_N function.
        """
        # TODO: need to add this correctly with the correct mass scales.
    

    def get_A(self):
        """
        Code equivalent of equation 2.7 in https://arxiv.org/abs/2209.01666.
        """
        # Note: replaces the A_func including fixing the second term which should be 2*G3phi rather than G3phi
        self.symfunc['A'] = (self.sym['M_Ks']**2.)*self.symfunc['Kx'] - 2*self.sym['M_G3s']*self.symfunc['G3phi'] + 2.*self.sym['X']*self.sym['M_G3s']*self.symfunc['G3phix']
        self.symfunc['A'] += 6*(self.sym['E']**2)*self.sym['phi_prime']*(self.sym['M_G3s']*self.symfunc['G3x'] + self.sym['X']*self.sym['M_G3s']*self.symfunc['G3xx']) 
        self.symfunc['A'] += (self.sym['E']**2.)*(self.sym['phi_prime']**2.)*((self.sym['M_Ks']**2.)*self.symfunc['Kxx'] - 2.*self.sym['M_G3s']*self.symfunc['G3phix'])


    def get_A_new(self):
        """
        Code equivalent of equation 2.7 in https://arxiv.org/abs/2209.01666.
        """
        self.symfunc['A_new'] = self.sym['M_Kp2']*self.symfunc['Kx']
        self.symfunc['A_new'] -= 2*self.sym['M_G3p']*self.sym['M_sp']*self.symfunc['G3phi']
        self.symfunc['A_new'] += (self.sym['E']**2)*6*self.sym['M_G3p']*self.sym['M_sp']*self.sym['phi_prime']*(self.symfunc['G3x']+self.sym['X']*self.symfunc['G3xx'])
        self.symfunc['A_new'] += (self.sym['phi_prime']**2)*(self.sym['M_Kp2']*self.symfunc['Kxx']-2*self.sym['M_G3p']*self.sym['M_sp']*self.symfunc['G3phix'])
    

    def get_B1(self):
        """
        Code equivalent of equation 2.9 in https://arxiv.org/abs/2209.01666.
        """
        self.symfunc['B1'] = 6*self.sym['X']*self.symfunc['G3x'] - 6*self.symfunc['G4phi']


    def get_B1_new(self):
        """
        Code equivalent of equation 2.9 in https://arxiv.org/abs/2209.01666.
        """
        self.symfunc['B1_new'] = 6*self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.symfunc['G3x'] - 6*self.sym['M_G4p2']*self.symfunc['G4phi']


    def get_B2(self):
        """
        Code equivalent of equation 2.10 in https://arxiv.org/abs/2209.01666.
        """
        self.symfunc['B2'] = 3*self.sym['phi_prime']*(self.symfunc['Kx'] - 2*self.symfunc['G3phi'] + 2*self.sym['X']*self.symfunc['G3phix'])
        self.symfunc['B2'] += (self.sym['phi_prime']**2)*(self.symfunc['Kxphi'] - 2*self.symfunc['G3phiphi'])
        self.symfunc['B2'] -= self.symfunc['Kphi']/(self.sym['E']**2) + 12*self.symfunc['G4phi']
        self.symfunc['B2'] += 18*self.sym['X']*self.symfunc['G3x'] + 2*self.sym['X']*self.symfunc['G3phiphi']/(self.sym['E']**2)
    

    def get_B2_new(self):
        """
        Code equivalent of equation 2.10 in https://arxiv.org/abs/2209.01666.
        """
        self.symfunc['B2_new'] = 3*self.sym['phi_prime']*(self.sym['M_Kp2']*self.symfunc['Kx']-2*self.sym['M_G3p']*self.sym['M_sp']*self.symfunc['G3phi'] + 2*self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.symfunc['G3xphi'])
        self.symfunc['B2_new'] += (self.sym['phi_prime']**2)*(self.sym['M_Kp2']*self.symfunc['Kxphi'] - 2*self.sym['M_G3p']*self.sym['M_sp']*self.symfunc['G3phiphi'])
        self.symfunc['B2_new'] -= self.sym['M_Kp2']*self.symfunc['Kphi']/(self.sym['E']**2) 
        self.symfunc['B2_new'] -= 12*self.sym['M_G3p']*self.symfunc['G4phi'] 
        self.symfunc['B2_new'] += 18*self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.symfunc['G3x'] 
        self.symfunc['B2_new'] += 2*self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*self.symfunc['G3phiphi']/(self.sym['E']**2)


    def get_C1_new(self):
        """
        TODO: need to add notes here
        """
        self.symfunc['C1_new'] = self.sym['M_Kp2']*self.symfunc['K'] 
        self.symfunc['C1_new'] -= 2*self.sym['M_G3p']*self.sym['M_sp']*self.sym['X']*(self.symfunc['G3phi']+(self.sym['E']*self.sym['E_prime']*self.sym['phi_prime'] + self.sym['E']**2 * self.sym['phi_primeprime'])*self.symfunc['G3x'])
        self.symfunc['C1_new'] += 2*self.sym['M_G4p2']*(self.sym['E']*self.sym['E_prime']*self.sym['phi_prime'] + self.sym['E']**2 * self.sym['phi_primeprime'] + 2*self.sym['E']**2 * self.sym['phi_prime'])*self.symfunc['G4phi'] 
        self.symfunc['C1_new'] += 4*self.sym['M_G4p2']*self.sym['X']*self.symfunc['G4phiphi']
        

    def get_C2_new(self):
        """
        TODO: need to add notes here
        """
        pressure_allnophi = (1./3.)*self.sym['Omega_g'] + self.sym['w_n1']*self.sym['Omega_n1'] + self.sym['w_n2']*self.sym['Omega_n2'] + self.sym['w_n3']*self.sym['Omega_n3'] + self.sym['w_l']*self.sym['Omega_l']
        self.symfunc['C2_new'] = 3*(self.sym['E']**2)*pressure_allnophi + 6*self.sym['M_G4p2']*(self.sym['E']**2)*self.symfunc['G4']
    

    def get_C_new(self):
        """
        TODO: need to add notes here
        """
        self.symfunc['C_new'] = self.symfunc['C1_new'] + self.symfunc['C2_new']


    def get_EprimeE(self):
        """
        Returns the evolution of the dimensionless Hubble function E=H/H0, as a function of log(a), follows equation 2.5 
        https://arxiv.org/pdf/2209.01666.pdf.
        """
        # Note: replaces EprimeEODERHS function.
        M_G4s = 1./self.sym['M_sG4']
        A = (self.sym['M_Ks']**2.)*self.symfunc['Kx'] - 2*self.sym['M_G3s']*self.symfunc['G3phi'] + 2.*self.sym['X']*self.sym['M_G3s']*self.symfunc['G3phix']
        A += 6*(self.sym['E']**2)*self.sym['phi_prime']*(self.sym['M_G3s']*self.symfunc['G3x'] + self.sym['X']*self.sym['M_G3s']*self.symfunc['G3xx']) 
        A += (self.sym['E']**2.)*(self.sym['phi_prime']**2.)*((self.sym['M_Ks']**2.)*self.symfunc['Kxx'] - 2.*self.sym['M_G3s']*self.symfunc['G3phix'])

        B1 = 6.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3x'] - 6.*(M_G4s**2)*self.symfunc['G4phi']

        B2 = 3.*self.sym['phi_prime']*((self.sym['M_Ks']**2.)*self.symfunc['Kx'] - 2.*self.sym['M_G3s']*self.symfunc['G3phi'] + 2.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3phix'])
        B2 += (self.sym['phi_prime']**2.)*((self.sym['M_Ks']**2.)*self.symfunc['Kxphi'] - 2.*self.sym['M_G3s']*self.symfunc['G3phiphi']) - ((self.sym['M_Ks']**2.)/(self.sym['E']**2.))*self.symfunc['Kphi']
        B2 += -12.*(M_G4s**2.)*self.symfunc['G4phi'] + 18.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3x'] + 2.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3phiphi']*(1./self.sym['E']**2.)

        term1 = 1.+ (1./2.)*self.sym['M_G3G4']*self.sym['M_sG4']*self.sym['X']*self.symfunc['G3x']*(1./self.symfunc['G4'])*(B1/A) - (self.symfunc['G4phi']*B1)/(2*self.symfunc['G4']*A)

        term21 = (self.sym['M_KG4']**2.)*self.symfunc['K']*(1./(self.sym['E']**2.)) - 2.*self.sym['M_sG4']*self.sym['M_G3G4']*(1./(self.sym['E']**2.))*self.sym['X']*self.symfunc['G3phi']
        term22 = 4.*self.symfunc['G4phi']*self.sym['phi_prime'] + 4.*self.sym['X']*self.symfunc['G4phiphi']*(1./(self.sym['E']**2.))
        term2 = (-1./(4.*self.symfunc['G4']))*(term21 + term22)

        term3 = (-1./2.)*((3*((1./3.)*self.sym['Omega_g'] + self.sym['w_n1']*self.sym['Omega_n1'] + self.sym['w_n2']*self.sym['Omega_n2'] + self.sym['w_n3']*self.sym['Omega_n3'] + self.sym['w_l']*self.sym['Omega_l'])/(2.*self.symfunc['G4']))*(self.sym['M_pG4']**2.) + 3.) 
        term4 = (-1./2.)*self.sym['M_G3G4']*self.sym['M_sG4']*self.sym['X']*self.symfunc['G3x']*(B2/(self.symfunc['G4']*A)) + (self.symfunc['G4phi']*B2)/(2*self.symfunc['G4']*A)
        RHS = term2 + term3 + term4

        self.symfunc['E_prime/E'] =  sym.simplify((RHS)/sym.simplify(term1))

    
    def get_Eprime_ODE_new(self):
        """
        Returns the Eprime ODE.
        """
        self.symfunc['E_prime_ODE_new'] = 4*self.sym['M_G4p2']*self.symfunc['G4']*self.sym['E']*self.sym['E_prime'] + self.symfunc['C_new']
    

    def get_phi_primeprime(self):
        """
        Return the phi_primeprime function, equation 2.6 in https://arxiv.org/abs/2209.01666
        """
        # Note: replaces the phi_primeprimeODERHS function
        M_G4s = 1./self.sym['M_sG4']
        A = (self.sym['M_Ks']**2.)*self.symfunc['Kx'] - self.sym['M_G3s']*self.symfunc['G3phi'] + 2.*self.sym['X']*self.sym['M_G3s']*self.symfunc['G3phix']
        A += 6*(self.sym['E']**2)*self.sym['phi_prime']*(self.sym['M_G3s']*self.symfunc['G3x'] + self.sym['X']*self.sym['M_G3s']*self.symfunc['G3xx']) + (self.sym['E']**2.)*(self.sym['phi_prime']**2.)*((self.sym['M_Ks']**2.)*self.symfunc['Kxx'] - 2.*self.sym['M_G3s']*self.symfunc['G3phix'])
        B1 = 6.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3x'] - 6.*(M_G4s**2)*self.symfunc['G4phi']
        B2 = 3.*self.sym['phi_prime']*((self.sym['M_Ks']**2.)*self.symfunc['Kx'] - 2.*self.sym['M_G3s']*self.symfunc['G3phi'] + 2.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3phix'])
        B2 += (self.sym['phi_prime']**2.)*((self.sym['M_Ks']**2.)*self.symfunc['Kxphi'] - 2.*self.sym['M_G3s']*self.symfunc['G3phiphi']) - ((self.sym['M_Ks']**2.)/(self.sym['E']**2.))*self.symfunc['Kphi']
        B2 += -12.*(M_G4s**2.)*self.symfunc['G4phi'] + 18.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3x'] + 2.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3phiphi']*(1./(self.sym['E']**2.))
        B =  B1*(self.sym['E_prime']/self.sym['E'])+B2
        # self.symfunc['phi_primeprime'] = -1.*((B/A) + (self.sym['E_prime']/self.sym['E'])*self.sym['phi_prime'])
        self.symfunc['phi_primeprime'] = -1.*((B/A) + (self.symfunc['E_prime/E'])*self.sym['phi_prime'])


    def get_phi_primeprime_ODE_new(self):
        """
        Returns the phi primeprime ODE.
        """
        self.symfunc['phi_primeprime_ODE_new'] = self.symfunc['A_new']*(self.sym['E']*self.sym['phi_primeprime'] + self.sym['E_prime']*self.sym['phi_prime']) + self.symfunc['B1_new']*self.sym['E_prime'] + self.symfunc['B2_new']*self.sym['E']


    def get_fried_closure(self):
        """
        Returns the RHS of the Friedmann closure relation -> equation 2.3 in https://arxiv.org/abs/2209.01666, 
        should be equal to 1. 
        """
        # Note: replaces the fried_closure function
        self.get_Omega_phi()
        self.symfunc['fried_closure'] = self.symfunc['Omega_phi'] + self.sym['Omega_b'] + self.sym['Omega_c'] + self.sym['Omega_g'] + self.sym['Omega_n1'] + self.sym['Omega_n2'] + self.sym['Omega_n3'] + self.sym['Omega_l'] - 1.
        # equation A.3 in https://arxiv.org/abs/2209.01666, 
        # TODO: duplication of a step done later in construct model, remove these lines.
        # Xreal = (1./2.)*(self.sym['E']**2.)*(self.sym['phi_prime']**2.)
        # self.symfunc['fried_closure'] = self.symfunc['fried_closure'].subs(self.sym['X'], Xreal)


    def get_theta(self):
        """
        Code equivalent of equation 2.15 in https://arxiv.org/abs/2209.01666
        """
        # Note: replaces theta function
        # All terms seemed to be multiplied by E, not so in eq. 2.15, why?
        term1 = self.sym['M_sG4']*self.sym['M_G3G4']*self.sym['E']*self.sym['phi_prime']*self.sym['X']*self.symfunc['G3x']/self.sym['M_pG4']/self.sym['M_pG4']
        term2 = 2.*self.sym['E']*self.symfunc['G4']/self.sym['M_pG4']/self.sym['M_pG4']
        term3 = self.sym['E']*self.sym['phi_prime']*self.symfunc['G4phi']/self.sym['M_pG4']/self.sym['M_pG4']
        self.symfunc['theta'] = -1.*term1 + term2 + term3
        Xreal = 0.5*(self.sym['E']**2.)*self.sym['phi_prime']**2.
        self.symfunc['theta'] = self.symfunc['theta'].subs(self.sym['X'], Xreal)
        

    def get_calE(self):
        """
        Code equivalent of equation 5 in https://arxiv.org/abs/1111.6749.
        """
        # Note: replaces calE function
        term1 = 2.*self.sym['M_KG4']*self.sym['M_KG4']*self.sym['X']*self.symfunc['Kx']/self.sym['M_pG4']/self.sym['M_pG4']
        term2 = self.sym['M_KG4']*self.sym['M_KG4']*self.symfunc['K']/self.sym['M_pG4']/self.sym['M_pG4']
        term3 = 6.*self.sym['M_sG4']*self.sym['M_G3G4']*self.sym['E']*self.sym['E']*self.sym['X']*self.symfunc['G3x']*self.sym['phi_prime']/self.sym['M_pG4']/self.sym['M_pG4']
        term4 = 2.*self.sym['M_sG4']*self.sym['M_G3G4']*self.sym['X']*self.symfunc['G3phi']/self.sym['M_pG4']/self.sym['M_pG4']
        term5 = 6.*self.sym['E']*self.sym['E']*self.symfunc['G4']/self.sym['M_pG4']/self.sym['M_pG4']
        term6 = 6.*self.sym['E']*self.sym['E']*self.symfunc['G4phi']*self.sym['phi_prime']/self.sym['M_pG4']/self.sym['M_pG4']
        self.symfunc['calE'] = term1 - term2 + term3 - term4 - term5 - term6


    def get_calP(self):
        """
        Code equivalent of equation 6 in https://arxiv.org/abs/1111.6749
        """
        # Note: replaces calP function
        term1 = self.sym['M_KG4']*self.sym['M_KG4']*self.symfunc['K']/self.sym['M_pG4']/self.sym['M_pG4']
        term2coeff = 2.*self.sym['M_sG4']*self.sym['X']/self.sym['M_pG4']
        term2bracket = self.sym['M_G3s']*self.symfunc['G3phi'] + self.sym['E']*self.sym['M_G3s']*self.symfunc['G3x']*(self.sym['E_prime']*self.sym['phi_prime'] + self.sym['E']*self.sym['phi_primeprime'])
        term2 = term2coeff*term2bracket
        term3 = 2.*self.symfunc['G4']*(3.*self.sym['E']*self.sym['E'] + 2*self.sym['E']*self.sym['E_prime'])/self.sym['M_pG4']/self.sym['M_pG4']
        term4coeff = 2*self.symfunc['G4phi']/self.sym['M_pG4']/self.sym['M_pG4']
        term4bracket = self.sym['E']*(self.sym['E_prime']*self.sym['phi_prime'] + self.sym['E']*self.sym['phi_primeprime']) + 2*self.sym['E']*self.sym['E']*self.sym['phi_prime']
        term4 = term4coeff*term4bracket
        self.symfunc['calP'] = term1 - term2 + term3 + term4


    def get_alpha0(self):
        """
        Code equivalent of equation 3.5 (see 2.16-2.19) in https://arxiv.org/abs/2209.01666.
        Given by alpha0 = A0/2G4.
        """
        # Note: replaces the alpha0 function

        self.get_calE()
        self.get_calP()
        self.get_theta()

        self.sym['dtheta_dE'] = sym.diff(self.symfunc['theta'], self.sym['E'])
        self.sym['dtheta_dphi_prime'] = sym.diff(self.symfunc['theta'], self.sym['phi_prime'])
        self.sym['thetaprime'] = self.sym['dtheta_dE']*self.sym['E_prime'] + self.sym['dtheta_dphi_prime']*self.sym['phi_primeprime']

        A0 = self.sym['thetaprime']/self.sym['E'] + self.symfunc['theta']/self.sym['E'] 
        A0 += -2.*self.symfunc['G4']/self.sym['M_pG4']/self.sym['M_pG4'] - 4.*self.symfunc['G4phi']*self.sym['phi_prime']/self.sym['M_pG4']/self.sym['M_pG4'] 
        A0 += -(self.symfunc['calE'] + self.symfunc['calP'])/(2.*self.sym['E']*self.sym['E'])
        self.symfunc['alpha0'] = self.sym['M_pG4']*self.sym['M_pG4']*A0/2./self.symfunc['G4']
        

    def get_alpha1(self):
        """
        Code equivalent of equation 3.5 (see 2.16-2.19) in https://arxiv.org/abs/2209.01666.
        Given by alpha1 = A1/2G4.
        """
        # Note: replaces the alpha1 function
        A1 = 2.*self.symfunc['G4phi']*self.sym['phi_prime']/self.sym['M_pG4']/self.sym['M_pG4']
        self.symfunc['alpha1'] = self.sym['M_pG4']*self.sym['M_pG4']*A1/2./self.symfunc['G4']


    def get_alpha2(self):
        """
        Code equivalent of equation 3.5 (see 2.16-2.19) in https://arxiv.org/abs/2209.01666.
        Given by alpha2 = A2/2G4.
        """
        # Note: replaces the alpha2 function
        self.get_theta()
        A2 = 2.*self.symfunc['G4']/self.sym['M_pG4']/self.sym['M_pG4'] - self.symfunc['theta']/self.sym['E']
        self.symfunc['alpha2'] = self.sym['M_pG4']*self.sym['M_pG4']*A2/2./self.symfunc['G4']


    def get_beta0(self):
        """
        Code equivalent of equation 3.5 (see 2.16-2.19) in https://arxiv.org/abs/2209.01666.
        Given by beta0 = B0/2G4.
        """
        # Note: replaces the beta0 function
        B0 = self.sym['M_sG4']*self.sym['M_G3G4']*self.sym['X']*self.symfunc['G3x']*self.sym['phi_prime']/self.sym['M_pG4']/self.sym['M_pG4']
        self.symfunc['beta0'] = self.sym['M_pG4']*self.sym['M_pG4']*B0/2./self.symfunc['G4']
    

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
    
    
    # Symbolic definitions for Bellini & Sawicki alphas

    def get_M_star_sq(self):
        """
        Code equivalent of equation A.6 in https://iopscience.iop.org/article/10.1088/1475-7516/2014/07/050 within reduced Horndeski class.
        """
        self.symfunc['M_star_sq'] = 2*self.symfunc['G4']*self.sym['M_p']**2
    

    def get_alpha_M(self):
        """
        Code equivalent of equation A.7 in https://iopscience.iop.org/article/10.1088/1475-7516/2014/07/050 within reduced Horndeski class.
        """
        self.symfunc['alpha_M'] = (2*self.sym['phi_prime']*self.symfunc['G4phi']*self.sym['M_p']**2)/self.symfunc['M_star_sq']
    

    def get_alpha_B(self):
        """
        Code equivalent of equation A.9 in https://iopscience.iop.org/article/10.1088/1475-7516/2014/07/050 within reduced Horndeski class.
        """
        self.symfunc['alpha_B'] = 2*self.sym['phi_prime']*(self.sym['M_s']*self.sym['M_G3']*self.sym['X']*self.symfunc['G3x'] - self.symfunc['G4phi']*self.sym['M_G4'])/self.symfunc['M_star_sq']


    def get_alpha_K(self):
        """
        Code equivalent of equation A.8 in https://iopscience.iop.org/article/10.1088/1475-7516/2014/07/050 within reduced Horndeski class.
        """
        term11 = (self.sym['M_K']**2)*(self.symfunc['Kx'] + 2*self.sym['X']*self.symfunc['Kxx'])
        term12 = 2*self.sym['M_s']*self.sym['M_G3']*(self.symfunc['G3phi'] + self.sym['X']*self.symfunc['G3phix'])
        term1 = (term11 - term12)*2*self.sym['X']/(self.symfunc['M_star_sq']*self.sym['E']**2)
        term2 = 12*self.sym['M_s']*self.sym['M_G3']*self.sym['phi_prime']*self.sym['X']*(self.symfunc['G3x'] + self.sym['X']*self.symfunc['G3xx'])/self.symfunc['M_star_sq']
        self.symfunc['alpha_K'] = term1 + term2
    
    # Symbolic definitions for scalar field density and pressure

    def get_rho_phi(self):
        """
        Code equivalent of equation A.1 in https://iopscience.iop.org/article/10.1088/1475-7516/2014/07/050 within reduced Horndeski class.
        Referred to as Tilde epsilon in internal Hi-COLA notes.
        """
        term1 = (self.sym['M_K']**2)*self.symfunc['Kx'] - self.sym['M_s']*self.sym['M_G3']*self.symfunc['G3phi']
        term2 = self.sym['M_s']*self.sym['M_G3']*self.symfunc['G3x']*self.sym['X'] - self.symfunc['G4phi']*self.sym['M_G4']**2
        # TODO: make this independent of H0, for easier H0 scaling
        self.symfunc['rho_phi'] = (self.sym['H0']**2)*(2*self.sym['X']*term1 + 6*(self.sym['E']**2)*self.sym['phi_prime']*term2 - self.sym['M_K']**2*self.symfunc['K'])/self.symfunc['M_star_sq']


    def get_P_phi(self):
        """
        Code equivalent of equation A.2 in https://iopscence.iop.org/article/10.1088/1475-7516/2014/07/050 within reduced Horndeski class.
        Referred to as Tilde P in internal Hi-COLA notes.
        """
        term1 = (self.sym['M_K']**2)*self.symfunc['K'] + 4*(self.sym['E']**2)*self.symfunc['G4phi']*self.sym['phi_prime']*(self.sym['M_G4']**2)
        term2 = self.sym['M_s']*self.sym['M_G3']*self.symfunc['G3phi'] - 2*self.symfunc['G4phiphi']*(self.sym['M_G4']**2)
        term31 = 2*self.sym['phi_prime']*(self.sym['M_s']*self.sym['M_G3']*self.sym['X']*self.symfunc['G3x'] - self.symfunc['G4phi']*(self.sym['M_G4']**2))
        term3 = term31*self.sym['E']*(self.sym['E_prime']*self.sym['phi_prime'] + self.sym['E']*self.sym['phi_primeprime'])/self.sym['phi_prime']
        # TODO: make this independent of H0, for easier H0 scaling
        self.symfunc['P_phi'] = (self.sym['H0']**2)*(term1 - 2*self.sym['X']*term2 - term3)/self.symfunc['M_star_sq']

    # Stability and sound speed related functions

    def get_Q_s(self):
        """
        Code equivalent of the LHS of the first inequality in 3.13 in https://iopscience.iop.org/article/10.1088/1475-7516/2014/07/050 within reduced Horndeski class.
        """
        self.symfunc['D'] = self.symfunc['alpha_K'] + (3/2)*(self.symfunc['alpha_B']**2)
        self.symfunc['Q_s'] = 2*self.symfunc['M_star_sq']*self.symfunc['D']/((2-self.symfunc['alpha_B'])**2) 
    

    def get_c_s_sq(self):
        """
        Code equivalent of the LHS of the second inequality in 3.13 in https://iopscience.iop.org/article/10.1088/1475-7516/2014/07/050 within reduced Horndeski class.
        """
        self.sym['alpha_B_prime'] = sym.symbols('alpha_{B}prime')
        term1 = self.symfunc['E_prime/E'] - self.symfunc['alpha_B']/2 - self.symfunc['alpha_M']
        term2 = 3*(self.sym['Omega_b'] + self.sym['Omega_c'] + (1.+1/3)*self.sym['Omega_g'] + (1. + self.sym['w_n1'])*self.sym['Omega_n1'] + (1. + self.sym['w_n2'])*self.sym['Omega_n2'] + (1. + self.sym['w_n3'])*self.sym['Omega_n3'] + (1. + self.sym['w_l'])*self.sym['Omega_l'])
        self.symfunc['c_s_sq_D'] = -((2 - self.symfunc['alpha_B'])*term1 - self.sym['alpha_B_prime'] + term2/self.symfunc['M_star_sq'])
        self.symfunc['c_s_sq'] = self.symfunc['c_s_sq_D']/self.symfunc['D']

    # Set mass ratios

    def set_mass_ratios(self, M_p=1, M_sp=1, M_gp=1, M_Kp=1, M_G3p=1, M_G4p=1):
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
        self.params['mass_ratios'] = {
            'M_p': M_p, 'M_s': M_sp/M_p, 'M_g': M_gp/M_p, 'M_K': M_Kp/M_p, 'M_G3': M_G3p/M_p, 'M_G4': M_G4p/M_p,  
        }
    
    def set_mass_ratios_new(self, M_sp=1, M_Kp2=1, M_G3p=1, M_G4p2=1):
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
        self.params['mass_ratios_new'] = {
            'M_sp': M_sp, 'M_Kp2': M_Kp2, 'M_G3p': M_G3p, 'M_G4p2': M_G4p2
        }


    def get_absolute_mass_ratios(self):
        """
        Allows the users to assign specific values to mass ratios, those unchanged will be set to 1.
        """
        if self.params['mass_ratios'] is None:
            self.set_mass_ratios()
        self.params['mass_ratios']['M_pG4'] = self.params['mass_ratios']['M_p']/self.params['mass_ratios']['M_G4']
        self.params['mass_ratios']['M_KG4'] = self.params['mass_ratios']['M_K']/self.params['mass_ratios']['M_G4']
        self.params['mass_ratios']['M_G3s'] = self.params['mass_ratios']['M_G3']/self.params['mass_ratios']['M_s']
        self.params['mass_ratios']['M_sG4'] = self.params['mass_ratios']['M_s']/self.params['mass_ratios']['M_G4']
        self.params['mass_ratios']['M_G3G4'] = self.params['mass_ratios']['M_G3']/self.params['mass_ratios']['M_G4']
        self.params['mass_ratios']['M_Ks'] = self.params['mass_ratios']['M_K']/self.params['mass_ratios']['M_s']
        self.params['mass_ratios']['M_gp'] = self.params['mass_ratios']['M_g']/self.params['mass_ratios']['M_p']
        self.params['mass_ratios']['M_sp'] = self.params['mass_ratios']['M_s']/self.params['mass_ratios']['M_p']


    def get_absolute_mass_ratios_new(self):
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
        if self.params['mass_ratios_new'] is None:
            self.set_mass_ratios_new()
    

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
                print(' - substituting X = 0.5 * E^2 * phi_prime^2')
                print(
                    ' - substituting mass ratios: M_p=%0.2f, M_s=%0.2f, M_g=%0.2f, M_K=%0.2f, M_G3=%0.2f, M_G4=%0.2f' % (
                        self.params['mass_ratios']['M_p'], self.params['mass_ratios']['M_s'], self.params['mass_ratios']['M_g'], 
                        self.params['mass_ratios']['M_K'], self.params['mass_ratios']['M_G3'], self.params['mass_ratios']['M_G4']
                    )
                )

            Xreal = 0.5*(self.sym['E']**2.)*self.sym['phi_prime']**2.

            sub_dict = {
                self.sym['X']: Xreal,
                self.sym['M_p']: self.params['mass_ratios']['M_p'],
                self.sym['M_pG4']: self.params['mass_ratios']['M_pG4'],
                self.sym['M_KG4']: self.params['mass_ratios']['M_KG4'],
                self.sym['M_G3s']: self.params['mass_ratios']['M_G3s'],
                self.sym['M_sG4']: self.params['mass_ratios']['M_sG4'],
                self.sym['M_G3G4']: self.params['mass_ratios']['M_G3G4'],
                self.sym['M_Ks']: self.params['mass_ratios']['M_Ks'],
                self.sym['M_gp']: self.params['mass_ratios']['M_gp'],
                self.sym['M_sp']: self.params['mass_ratios']['M_sp'],
                self.sym['M_G4']: self.params['mass_ratios']['M_G4'],
                self.sym['M_K']: self.params['mass_ratios']['M_K'],
                self.sym['M_s']: self.params['mass_ratios']['M_s'],
                self.sym['M_g']: self.params['mass_ratios']['M_g'],
                self.sym['M_G3']: self.params['mass_ratios']['M_G3'],
            }

            if self.verbose:
                print(' - into symbolic functions...')

            self.get_G_G_4()
            G_G_4_GN = self.symfunc['G_G_4/G_N'].subs(sub_dict)

            self.get_EprimeE()
            EprimeE = self.symfunc['E_prime/E'].subs(sub_dict)

            self.get_phi_primeprime()
            phi_primeprime = self.symfunc['phi_primeprime'].subs(sub_dict)

            self.get_A()
            A = self.symfunc['A'].subs(sub_dict)

            self.get_B1()
            B1 = self.symfunc['B1'].subs(sub_dict)

            self.get_B2()
            B2 = self.symfunc['B2'].subs(sub_dict)

            self.get_Omega_phi()
            Omega_phi = self.symfunc['Omega_phi'].subs(sub_dict)

            self.get_fried_closure()
            fried_closure = self.symfunc['fried_closure'].subs(sub_dict)

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
            
            self.get_rho_phi()
            rho_phi = self.symfunc['rho_phi'].subs(sub_dict)

            self.get_P_phi()
            P_phi = self.symfunc['P_phi'].subs(sub_dict)

            self.get_Q_s()
            D = self.symfunc['D'].subs(sub_dict)
            Q_s = self.symfunc['Q_s'].subs(sub_dict)

            self.get_c_s_sq()
            c_s_sq_D = self.symfunc['c_s_sq_D'].subs(sub_dict)
            c_s_sq = self.symfunc['c_s_sq'].subs(sub_dict)
            
            # we will copy these substituted and simplified functions to the class symfunc dictionary, we avoided 
            # doing this before as some of these are re-called and redefined in the 'get' functions.
            self.symfunc['G_G_4/G_N'] = G_G_4_GN
            self.symfunc['E_prime/E'] = EprimeE
            self.symfunc['phi_primeprime'] = phi_primeprime
            self.symfunc['A'] = A
            self.symfunc['B1'] = B1
            self.symfunc['B2'] = B2
            self.symfunc['Omega_phi'] = Omega_phi
            self.symfunc['fried_closure'] = fried_closure
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
            self.symfunc['alpha_K'] = alpha_K
            self.symfunc['rho_phi'] = rho_phi
            self.symfunc['P_phi'] = P_phi
            self.symfunc['D'] = D
            self.symfunc['Q_s'] = Q_s
            self.symfunc['c_s_sq_D'] = c_s_sq_D
            self.symfunc['c_s_sq'] = c_s_sq

            # Lambdify functions

            if self.verbose:
                print(" - 'Lambdify'ing symbolic functions")

            # keep variables fixed to functions solved in the ODE + Horndeski variables
            variables = [
                self.sym['E'], self.sym['phi'], self.sym['phi_prime'], self.sym['Omega_g'], self.sym['Omega_b'], self.sym['Omega_c'], self.sym['Omega_l'],
                self.sym['Omega_n1'], self.sym['w_n1'], self.sym['Omega_n2'], self.sym['w_n2'], self.sym['Omega_n3'], self.sym['w_n3'], self.sym['w_l'], 
                *self.sym['K_G3_G4_syms'], self.sym['f_H']]
            
            # Need to figure out how to input neutrino EoS, at the moment this is an effective EoS not a single one as calculated before.

            # G_G_4/G_N function
            self.lambda_funcs['G_G_4/G_N'] = sym.lambdify(variables, self.symfunc['G_G_4/G_N'], 'numpy')

            # functions as part of the ODE set of equations relating phi and E.
            self.lambda_funcs['B2_lambda'] = sym.lambdify(variables, self.symfunc['B2'], 'numpy')
            self.lambda_funcs['fried_closure_lambda'] = sym.lambdify(variables, self.symfunc['fried_closure'], 'numpy')
            self.lambda_funcs['E_prime/E_lambda'] = sym.lambdify(variables, self.symfunc['E_prime/E'], 'numpy')

            # Add E_prime to variables compute phi_primeprime
            variables = [self.sym['E_prime'], *variables]
            self.lambda_funcs['phi_primeprime_lambda'] = sym.lambdify(variables, self.symfunc['phi_primeprime'], 'numpy')
            
            # Derived quantities
            # Add phi_primeprime to variables to compute the rest of the derived quantities
            variables = [self.sym['phi_primeprime'], *variables]

            self.lambda_funcs['Omega_phi_lambda'] = sym.lambdify(variables, self.symfunc['Omega_phi'], 'numpy')
            self.lambda_funcs['A_lambda'] = sym.lambdify(variables, self.symfunc['A'], 'numpy')
            self.lambda_funcs['B1_lambda'] = sym.lambdify(variables, self.symfunc['B1'], 'numpy')
            self.lambda_funcs['B2_lambda'] = sym.lambdify(variables, self.symfunc['B2'], 'numpy')
            self.lambda_funcs['alpha0_lambda'] = sym.lambdify(variables, self.symfunc['alpha0'], 'numpy')
            self.lambda_funcs['alpha1_lambda'] = sym.lambdify(variables, self.symfunc['alpha1'], 'numpy')
            self.lambda_funcs['alpha2_lambda'] = sym.lambdify(variables, self.symfunc['alpha2'], 'numpy')
            self.lambda_funcs['beta0_lambda'] = sym.lambdify(variables, self.symfunc['beta0'], 'numpy')
            self.lambda_funcs['calB_lambda'] = sym.lambdify(variables, self.symfunc['calB'], 'numpy')
            self.lambda_funcs['calC_lambda'] = sym.lambdify(variables, self.symfunc['calC'], 'numpy')
            self.lambda_funcs['beta_lambda'] = sym.lambdify(variables, self.symfunc['beta'], 'numpy')
            
            self.lambda_funcs['M_star_sq'] = sym.lambdify(variables, self.symfunc['M_star_sq'], 'numpy')
            self.lambda_funcs['alpha_M'] = sym.lambdify(variables, self.symfunc['alpha_M'], 'numpy')
            self.lambda_funcs['alpha_B'] = sym.lambdify(variables, self.symfunc['alpha_B'], 'numpy')
            self.lambda_funcs['alpha_K'] = sym.lambdify(variables, self.symfunc['alpha_K'], 'numpy')
            self.lambda_funcs['rho_phi'] = sym.lambdify([self.sym['H0'], *variables], self.symfunc['rho_phi'], 'numpy')
            self.lambda_funcs['P_phi'] = sym.lambdify([self.sym['H0'], *variables], self.symfunc['P_phi'], 'numpy')
            self.lambda_funcs['D'] = sym.lambdify(variables, self.symfunc['D'], 'numpy')
            self.lambda_funcs['Q_s'] = sym.lambdify(variables, self.symfunc['Q_s'], 'numpy')
            self.lambda_funcs['c_s_sq_D'] = sym.lambdify([self.sym['alpha_B_prime'], *variables], self.symfunc['c_s_sq_D'], 'numpy')
            self.lambda_funcs['c_s_sq'] = sym.lambdify([self.sym['alpha_B_prime'], *variables], self.symfunc['c_s_sq'], 'numpy')

            if self.verbose:
                print(' - Done!')


    def construct_model_new(self):
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

            self.get_absolute_mass_ratios_new()

            if self.verbose:
                print(' - substituting X = 0.5 * E^2 * phi_prime^2')
                print(
                    ' - substituting mass ratios: M_sp=%0.2f, M_Kp2=%0.2f, M_G3p=%0.2f, M_G4p2=%0.2f' % (
                        self.params['mass_ratios_new']['M_sp'], self.params['mass_ratios_new']['M_Kp2'], 
                        self.params['mass_ratios_new']['M_G3p'], self.params['mass_ratios_new']['M_G4p2'],
                    )
                )

            Xreal = 0.5*(self.sym['E']**2.)*self.sym['phi_prime']**2.

            sub_dict = {
                self.sym['X']: Xreal,
                self.sym['M_sp']: self.params['mass_ratios_new']['M_sp'],
                self.sym['M_Kp2']: self.params['mass_ratios_new']['M_Kp2'],
                self.sym['M_G3p']: self.params['mass_ratios_new']['M_G3p'],
                self.sym['M_G4p2']: self.params['mass_ratios_new']['M_G4p2'],
            }

            if self.verbose:
                print(' - into symbolic functions...')

            # self.get_G_G_4()
            # G_G_4_GN = self.symfunc['G_G_4/G_N'].subs(sub_dict)
            
            # self.get_A()
            # A = self.symfunc['A'].subs(sub_dict)

            # self.get_B1()
            # B1 = self.symfunc['B1'].subs(sub_dict)

            # self.get_B2()
            # B2 = self.symfunc['B2'].subs(sub_dict)

            self.get_Omega_phi_new()
            Omega_phi = self.symfunc['Omega_phi_new'].subs(sub_dict)

            self.get_fried_closure_new()
            fried_closure = self.symfunc['fried_closure_new'].subs(sub_dict)

            self.get_A_new()
            A = self.symfunc['A_new'].subs(sub_dict)

            self.get_B1_new()
            B1 = self.symfunc['B1_new'].subs(sub_dict)

            self.get_B2_new()
            B2 = self.symfunc['B2_new'].subs(sub_dict)

            self.get_C1_new()
            self.get_C2_new()
            self.get_C_new()
            C = self.symfunc['C_new'].subs(sub_dict)

            self.get_Eprime_ODE_new()
            E_prime_ODE = self.symfunc['E_prime_ODE_new'].subs(sub_dict)

            self.get_phi_primeprime_ODE_new()
            phi_primeprime_ODE = self.symfunc['phi_primeprime_ODE_new'].subs(sub_dict)

            # self.get_alpha0()
            # alpha0 = self.symfunc['alpha0'].subs(sub_dict)

            # self.get_alpha1()
            # alpha1 = self.symfunc['alpha1'].subs(sub_dict)

            # self.get_alpha2()
            # alpha2 = self.symfunc['alpha2'].subs(sub_dict)

            # self.get_beta0()
            # beta0 = self.symfunc['beta0'].subs(sub_dict)

            # self.get_calB()
            # calB = self.symfunc['calB'].subs(sub_dict)

            # self.get_calC()
            # calC = self.symfunc['calC'].subs(sub_dict)

            # self.get_beta()
            # beta = self.symfunc['beta'].subs(sub_dict)

            # self.get_M_star_sq()
            # M_star_sq = self.symfunc['M_star_sq'].subs(sub_dict)
            
            # self.get_alpha_M()
            # alpha_M = self.symfunc['alpha_M'].subs(sub_dict)

            # self.get_alpha_B()
            # alpha_B = self.symfunc['alpha_B'].subs(sub_dict)

            # self.get_alpha_K()
            # alpha_K = self.symfunc['alpha_K'].subs(sub_dict)
            
            # self.get_rho_phi()
            # rho_phi = self.symfunc['rho_phi'].subs(sub_dict)

            # self.get_P_phi()
            # P_phi = self.symfunc['P_phi'].subs(sub_dict)

            # self.get_Q_s()
            # D = self.symfunc['D'].subs(sub_dict)
            # Q_s = self.symfunc['Q_s'].subs(sub_dict)

            # self.get_c_s_sq()
            # c_s_sq_D = self.symfunc['c_s_sq_D'].subs(sub_dict)
            # c_s_sq = self.symfunc['c_s_sq'].subs(sub_dict)
            
            # we will copy these substituted and simplified functions to the class symfunc dictionary, we avoided 
            # doing this before as some of these are re-called and redefined in the 'get' functions.
            self.symfunc['E_prime_ODE'] = E_prime_ODE
            self.symfunc['phi_primeprime_ODE'] = phi_primeprime_ODE
            # self.symfunc['G_G_4/G_N'] = G_G_4_GN
            # self.symfunc['E_prime/E'] = EprimeE
            # self.symfunc['phi_primeprime'] = phi_primeprime
            # self.symfunc['A'] = A
            # self.symfunc['B1'] = B1
            # self.symfunc['B2'] = B2
            self.symfunc['Omega_phi'] = Omega_phi
            self.symfunc['fried_closure'] = fried_closure
            # self.symfunc['alpha0'] = alpha0
            # self.symfunc['alpha1'] = alpha1
            # self.symfunc['alpha2'] = alpha2
            # self.symfunc['beta0'] = beta0
            # self.symfunc['calB'] = calB
            # self.symfunc['calC'] = calC
            # self.symfunc['beta'] = beta
            # self.symfunc['M_star_sq'] = M_star_sq
            # self.symfunc['alpha_M'] = alpha_M
            # self.symfunc['alpha_B'] = alpha_B
            # self.symfunc['alpha_K'] = alpha_K
            # self.symfunc['rho_phi'] = rho_phi
            # self.symfunc['P_phi'] = P_phi
            # self.symfunc['D'] = D
            # self.symfunc['Q_s'] = Q_s
            # self.symfunc['c_s_sq_D'] = c_s_sq_D
            # self.symfunc['c_s_sq'] = c_s_sq

            # Lambdify functions

            if self.verbose:
                print(" - 'Lambdify'ing symbolic functions")

            # keep variables fixed to functions solved in the ODE + Horndeski variables
            variables = [
                self.sym['E'], self.sym['E_prime'], self.sym['phi'], self.sym['phi_prime'], self.sym['phi_primeprime'], 
                self.sym['Omega_g'], self.sym['Omega_b'], self.sym['Omega_c'], self.sym['Omega_l'],
                self.sym['Omega_n1'], self.sym['w_n1'], self.sym['Omega_n2'], self.sym['w_n2'], self.sym['Omega_n3'], self.sym['w_n3'], self.sym['w_l'], 
                *self.sym['K_G3_G4_syms'], self.sym['f_H']]
            
            # Need to figure out how to input neutrino EoS, at the moment this is an effective EoS not a single one as calculated before.

            # G_G_4/G_N function
            # self.lambda_funcs['G_G_4/G_N'] = sym.lambdify(variables, self.symfunc['G_G_4/G_N'], 'numpy')

            # # functions as part of the ODE set of equations relating phi and E.
            self.lambda_funcs['E_prime_ODE'] = sym.lambdify(variables, self.symfunc['E_prime_ODE'], 'numpy')
            self.lambda_funcs['phi_primeprime_ODE'] = sym.lambdify(variables, self.symfunc['phi_primeprime_ODE'], 'numpy')

            # self.lambda_funcs['B2_lambda'] = sym.lambdify(variables, self.symfunc['B2'], 'numpy')
            self.lambda_funcs['fried_closure_lambda'] = sym.lambdify(variables, self.symfunc['fried_closure'], 'numpy')
            # self.lambda_funcs['E_prime/E_lambda'] = sym.lambdify(variables, self.symfunc['E_prime/E'], 'numpy')

            # # Add E_prime to variables compute phi_primeprime
            # variables = [self.sym['E_prime'], *variables]
            # self.lambda_funcs['phi_primeprime_lambda'] = sym.lambdify(variables, self.symfunc['phi_primeprime'], 'numpy')
            
            # # Derived quantities
            # # Add phi_primeprime to variables to compute the rest of the derived quantities
            # variables = [self.sym['phi_primeprime'], *variables]

            self.lambda_funcs['Omega_phi_lambda'] = sym.lambdify(variables, self.symfunc['Omega_phi_new'], 'numpy')
            # self.lambda_funcs['A_lambda'] = sym.lambdify(variables, self.symfunc['A'], 'numpy')
            # self.lambda_funcs['B1_lambda'] = sym.lambdify(variables, self.symfunc['B1'], 'numpy')
            # self.lambda_funcs['B2_lambda'] = sym.lambdify(variables, self.symfunc['B2'], 'numpy')
            # self.lambda_funcs['alpha0_lambda'] = sym.lambdify(variables, self.symfunc['alpha0'], 'numpy')
            # self.lambda_funcs['alpha1_lambda'] = sym.lambdify(variables, self.symfunc['alpha1'], 'numpy')
            # self.lambda_funcs['alpha2_lambda'] = sym.lambdify(variables, self.symfunc['alpha2'], 'numpy')
            # self.lambda_funcs['beta0_lambda'] = sym.lambdify(variables, self.symfunc['beta0'], 'numpy')
            # self.lambda_funcs['calB_lambda'] = sym.lambdify(variables, self.symfunc['calB'], 'numpy')
            # self.lambda_funcs['calC_lambda'] = sym.lambdify(variables, self.symfunc['calC'], 'numpy')
            # self.lambda_funcs['beta_lambda'] = sym.lambdify(variables, self.symfunc['beta'], 'numpy')
            
            # self.lambda_funcs['M_star_sq'] = sym.lambdify(variables, self.symfunc['M_star_sq'], 'numpy')
            # self.lambda_funcs['alpha_M'] = sym.lambdify(variables, self.symfunc['alpha_M'], 'numpy')
            # self.lambda_funcs['alpha_B'] = sym.lambdify(variables, self.symfunc['alpha_B'], 'numpy')
            # self.lambda_funcs['alpha_K'] = sym.lambdify(variables, self.symfunc['alpha_K'], 'numpy')
            # self.lambda_funcs['rho_phi'] = sym.lambdify([self.sym['H0'], *variables], self.symfunc['rho_phi'], 'numpy')
            # self.lambda_funcs['P_phi'] = sym.lambdify([self.sym['H0'], *variables], self.symfunc['P_phi'], 'numpy')
            # self.lambda_funcs['D'] = sym.lambdify(variables, self.symfunc['D'], 'numpy')
            # self.lambda_funcs['Q_s'] = sym.lambdify(variables, self.symfunc['Q_s'], 'numpy')
            # self.lambda_funcs['c_s_sq_D'] = sym.lambdify([self.sym['alpha_B_prime'], *variables], self.symfunc['c_s_sq_D'], 'numpy')
            # self.lambda_funcs['c_s_sq'] = sym.lambdify([self.sym['alpha_B_prime'], *variables], self.symfunc['c_s_sq'], 'numpy')

            if self.verbose:
                print(' - Done!')
    
    # Explore scaling symmetry

    def scaling_symmetry(self, K_G3_G4_sub):
        """
        Allows the user to directly see how scaling symmetries apply for a given model.
        """
        sub_dict = {
            self.sym['f_H']: 1,
            self.sym['X']: sym.symbols('C_{\phi}')*sym.symbols('C_{\phi}')*self.sym['X'],
            self.sym['phi']: sym.symbols('C_{\phi}')*sym.symbols('\phi'),
            self.sym['phi_prime']: sym.symbols('C_{\phi}')*sym.symbols("\phi^{'}"),
            self.sym['phi_primeprime']: sym.symbols('C_{\phi}')*sym.symbols("\phi^{''}")
        }
        for (i, var) in enumerate(self.sym['K_G3_G4_syms']):
            # sub_dict[var] = sym.symbols(K_G3_G4_sub[i])
            sub_dict[var] = K_G3_G4_sub[i]
        
        E_primeE_term1 = self.symfunc['K']/(self.sym['E']**2)
        E_primeE_term1 -= 2*self.sym['X']*(self.symfunc['G3phi']/(self.sym['E']**2) + self.symfunc['G3x']*(self.sym['phi_prime']*self.sym['E_prime']/self.sym['E'] + self.sym['phi_primeprime']))
        E_primeE_term1 += 2*self.symfunc['G4phi']*(self.sym['phi_prime']*(self.sym['E_prime']/self.sym['E']+2) + self.sym['phi_primeprime'])
        E_primeE_term1 += 4*self.sym['X']*self.symfunc['G4phiphi']/(self.sym['E']**2)
        E_primeE_term1 /= 4*self.symfunc['G4']
        E_primeE_term2 = - ((self.sym['Omega_g'] - 3*self.sym['Omega_l'] + 3*self.sym['w_n1']*self.sym['Omega_n1'] + 3*self.sym['w_n2']*self.sym['Omega_n2'] + 3*self.sym['w_n3']*self.sym['Omega_n3'])/(2*self.symfunc['G4']) + 3)/2
        E_primeE = E_primeE_term1 + E_primeE_term2
        E_primeE = E_primeE.subs(sub_dict)

        A_term1 = 6*self.sym['phi_prime']*(self.symfunc['G3x']+self.sym['X']*self.symfunc['G3xx'])
        A_term2 = (self.sym['phi_prime']**2)*(self.symfunc['Kxx'] - 2*self.symfunc['G3phix'])
        A = self.symfunc['Kx'] - self.symfunc['G3phi'] + 2*self.sym['X']*self.symfunc['G3phix'] + (self.sym['E']**2)*(A_term1 + A_term2)

        B1 = 6*self.sym['X']*self.symfunc['G3x'] - 6*self.symfunc['G4phi']
        B2_term1 = 3*self.sym['phi_prime']*(self.symfunc['Kx']-2*self.symfunc['G3phi']+2*self.sym['X']*self.symfunc['G3phix'])
        B2_term2 = (self.sym['phi_prime']**2)*(self.symfunc['Kxphi'] - 2*self.symfunc['G3phiphi'])
        B2_term3 = -self.symfunc['Kphi']/(self.sym['E']**2) - 12*self.symfunc['G4phi'] + 18*self.sym['X']*self.symfunc['G3x'] + 2*self.sym['X']*self.symfunc['G3phiphi']/(self.sym['E']**2)
        B2 = B2_term1 + B2_term2 + B2_term3
        B = B1*(self.sym['E_prime']/self.sym['E']) + B2

        C_phi_phi_primeprime = - B/A - (self.sym['E_prime']/self.sym['E'])*self.sym['phi_prime']
        C_phi_phi_primeprime = C_phi_phi_primeprime.subs(sub_dict)

        fried_closure = self.sym['Omega_g'] + self.sym['Omega_c']  + self.sym['Omega_b'] + self.sym['Omega_l'] + self.sym['Omega_n1'] + self.sym['Omega_n2'] + self.sym['Omega_n3']
        Omega_phi_term1 = (self.sym['Omega_g'] + self.sym['Omega_c'] + self.sym['Omega_b'] + self.sym['Omega_l'] + self.sym['Omega_n1'] + self.sym['Omega_n2'] + self.sym['Omega_n3'])*(1/(2*self.symfunc['G4']) - 1)
        Omega_phi_term2 = self.sym['X']*self.symfunc['Kx']/(self.sym['E']**2)
        Omega_phi_term2 -= self.symfunc['K']/(2*self.sym['E']**2)
        Omega_phi_term2 += 3*self.sym['X']*self.sym['phi_prime']*self.symfunc['G3x']
        Omega_phi_term2 -= self.sym['X']*self.symfunc['G3phi']/(self.sym['E']**2)
        Omega_phi_term2 -= 3*self.sym['phi_prime']*self.symfunc['G4phi']
        Omega_phi_term2 *= (1/(3*self.symfunc['G4']))
        Omega_phi = Omega_phi_term1 + Omega_phi_term2
        fried_closure += Omega_phi

        fried_closure = fried_closure.subs(sub_dict)

        return E_primeE, C_phi_phi_primeprime, fried_closure
    
    
    # Densities for standard components of the Universe.

    ## Photons are given by the CMB temperature

    def _get_C_gamma(self):
        """
        Retrieve C_gamma constant.
        """
        self.const['C_gamma'] =  8 * (np.pi**3) * self.const['G[eV]'] * (self.params['Tcmb0'] * self.const['kB[eV]'])**4
        self.const['C_gamma'] /= 45*(self.const['H0unit[eV]']**2)
    

    def get_Omega_g_E2(self, a, h):
        """
        Returns the E^2 Omega_gamma (photon) density.

        Parameters
        ----------
        a : float or array
            Scale factor.
        h : float
            Little h = H0 * 1e-2.
        """
        return self.const['C_gamma']/((h**2)*(a**4))
    

    def get_Omega_g_E2_prime(self, a, h):
        """
        Returns the E^2 Omega_gamma (photon) density derivative.

        Parameters
        ----------
        a : float or array
            Scale factor.
        h : float
            Little h = H0 * 1e-2.
        """
        return -4*self.get_Omega_g_E2(a, h)
    


    def get_Omega_g(self, a, h, E):
        """
        Returns the Omega_g (photon) fractional density.

        Parameters
        ----------
        a : float or array
            Scale factor.
        h : float
            Little h = H0 * 1e-2.
        E : float or array
            The normalised Hubble expansion rate.
        """
        return self.get_Omega_g_E2(a, h)/(E**2)


    def compute_Omega_g_prime(self, Omega_g, E, E_prime):
        """
        Computes Omega photon prime.

        Parameters
        ----------
        Omega_g : float or array
            Radiation density.
        E : float or array
            Normalised Hubble expansion.
        E_prime : float or array
            Derivative of the normalised Hubble expansion.
        
        Returns
        -------
        Omega_g_prime : float or array
            Omega radiation prime.
        """
        E_prime_E = E_prime/E
        Omega_g_prime = -Omega_g*(4. + 2.*E_prime_E)
        return Omega_g_prime

    
    ## Neutrinos, given by CMB temperature and Neutrino mass and Hierarchy.

    def _get_C_nu(self):
        """
        Retrieve C_nu constant.
        """
        # Conventionally the Neutrino temperature would be given by the CMB temperature via
        # self.params['Tnu0'] = ((4/11)**(1/3))*self.params['Tcmb0']
        # however, in practice there are small corrections due to QED, e-e+ annihilation which means the neutrino temperature
        # is not quite this. So instead we use a value that reproduces the widely used mv/93.14h^2 approximation used.
        self.const['C_nu'] = (8 * self.const['G[eV]'] * (self.params['Tnu0'] * self.const['kB[eV]'])**4)/(3.*(self.const['H0unit[eV]']**2)*np.pi)
    

    def _tabulate_IJy(self):
        """
        Construct tabulated and interpolation function for Iy and Iy_prime used to compute neutrino density evolution.
        """

        from scipy.integrate import quad
        from scipy.interpolate import interp1d

        self.neutrino_table = {}

        def get_Iy(x, y):
            f = (x**2)*np.sqrt(x**2 + y**2) * np.exp(-x) / (1 + np.exp(-x))
            return f

        def get_Iy_prime(x, y):
            f = ((y**2)*(x**2)/np.sqrt(x**2 + y**2)) * np.exp(-x) / (1 + np.exp(-x))
            return f
        
        def get_Jy(x, y):
            f = ((x**4)/np.sqrt(x**2 + y**2)) * np.exp(-x) / (1 + np.exp(-x))
            return f
        
        y = np.logspace(-1, 2, 100)
        ymin = y.min()
        ymax = y.max()

        Iy = np.array([quad(get_Iy, 0., np.inf, args=(_y))[0] for _y in y])

        Iy_prime = np.array([quad(get_Iy_prime, 0., np.inf, args=(_y))[0] for _y in y])

        Jy = np.array([quad(get_Jy, 0., np.inf, args=(_y))[0] for _y in y])

        self.neutrino_table['y'] = y
        self.neutrino_table['ymin'] = ymin
        self.neutrino_table['ymax'] = ymax
        self.neutrino_table['Iy'] = Iy
        self.neutrino_table['Iy_prime'] = Iy_prime
        self.neutrino_table['Jy'] = Jy
        self.neutrino_table['Iy_interp'] = interp1d(self.neutrino_table['y'], self.neutrino_table['Iy'], kind='cubic')
        self.neutrino_table['Iy_prime_interp'] = interp1d(self.neutrino_table['y'], self.neutrino_table['Iy_prime'], kind='cubic')
        self.neutrino_table['Jy_interp'] = interp1d(self.neutrino_table['y'], self.neutrino_table['Jy'], kind='cubic')

        def Iy_asymp_low(y):
            return 7*(np.pi**4)/120
        
        self.neutrino_table['Iy_interp_low'] = Iy_asymp_low

        def Iy_asymp_high(y):
            return 3*y*self.const['riemann_zeta3']/2
        
        self.neutrino_table['Iy_interp_high'] = Iy_asymp_high

        def Iy_prime_asymp_low(y):
            return (np.pi**2) * (y**2) / 12
        
        self.neutrino_table['Iy_prime_interp_low'] = Iy_prime_asymp_low

        def Iy_prime_asymp_high(y):
            return 3*y*self.const['riemann_zeta3']/2
        
        self.neutrino_table['Iy_prime_interp_high'] = Iy_prime_asymp_high

        def Jy_asymp_low(y):
            return 7*(np.pi**4)/120
        
        self.neutrino_table['Jy_interp_low'] = Jy_asymp_low

        def Jy_asymp_high(y):
            return 45*self.const['riemann_zeta5']/(2*y)
        
        self.neutrino_table['Jy_interp_high'] = Jy_asymp_high


    def get_Omega_nu_E2(self, a, h, mnu):
        """
        Returns the neutrino density.

        Parameters
        ----------
        a : float or array
            Scale factor.
        h : float
            Little h = H0 * 1e-2.
        mnu : float
            Sets the mass of a single neutrino species.
        
        Returns
        -------
        Omega_nu_E2 : float or array
            Neutrino density.
        """
        if self.neutrino_table is None:
            self._tabulate_IJy()
        y = a * mnu / (self.params['Tnu0']*self.const['kB[eV]'])
        if np.isscalar(y):
            if y <= self.neutrino_table['ymin']:
                Omega_nu_E2 = self.neutrino_table['Iy_interp_low'](y)
            elif y >= self.neutrino_table['ymax']:
                Omega_nu_E2 = self.neutrino_table['Iy_interp_high'](y)
            else:
                Omega_nu_E2 = self.neutrino_table['Iy_interp'](y)
        else:
            Omega_nu_E2 = np.zeros(len(y))
            cond = np.where((y <= self.neutrino_table['ymin']))[0]
            Omega_nu_E2[cond] = self.neutrino_table['Iy_interp_low'](y[cond])
            cond = np.where((y >= self.neutrino_table['ymax']))[0]
            Omega_nu_E2[cond] = self.neutrino_table['Iy_interp_high'](y[cond])
            cond = np.where((y > self.neutrino_table['ymin']) & (y < self.neutrino_table['ymax']))[0]
            Omega_nu_E2[cond] = self.neutrino_table['Iy_interp'](y[cond])
        Omega_nu_E2 *= self.const['C_nu']/((h**2) * (a**4))
        return Omega_nu_E2
    
    
    def get_Omega_nu_E2_prime(self, a, h, mnu):
        """
        Returns the neutrino density prime.

        Parameters
        ----------
        a : float or array
            Scale factor.
        h : float
            Little h = H0 * 1e-2.
        mnu : float
            Sets the mass of a single neutrino species.
        
        Returns
        -------
        Omega_nu_E2_prime : float or array
            Neutrino density.
        """
        Omega_nu_E2_prime = -4*self.get_Omega_nu_E2(a, h, mnu)
        y = a * mnu / (self.params['Tnu0']*self.const['kB[eV]'])
        if np.isscalar(y):
            if y <= self.neutrino_table['ymin']:
                _Omega_nu_E2_prime = self.neutrino_table['Iy_prime_interp_low'](y)
            elif y >= self.neutrino_table['ymax']:
                _Omega_nu_E2_prime = self.neutrino_table['Iy_prime_interp_high'](y)
            else:
                _Omega_nu_E2_prime = self.neutrino_table['Iy_prime_interp'](y)
        else:
            _Omega_nu_E2_prime = np.zeros(len(y))
            cond = np.where((y <= self.neutrino_table['ymin']))[0]
            _Omega_nu_E2_prime[cond] = self.neutrino_table['Iy_prime_interp_low'](y[cond])
            cond = np.where((y >= self.neutrino_table['ymax']))[0]
            _Omega_nu_E2_prime[cond] = self.neutrino_table['Iy_prime_interp_high'](y[cond])
            cond = np.where((y > self.neutrino_table['ymin']) & (y < self.neutrino_table['ymax']))[0]
            _Omega_nu_E2_prime[cond] = self.neutrino_table['Iy_prime_interp'](y[cond])
        _Omega_nu_E2_prime *= self.const['C_nu']/((h**2) * (a**4))
        Omega_nu_E2_prime += _Omega_nu_E2_prime
        return Omega_nu_E2_prime
    
    
    def get_Omega_nu(self, a, h, mnu, E):
        """
        Returns the Omega_nu fractional density for a specific neutrinos.

        Parameters
        ----------
        a : float or array
            Scale factor.
        h : float
            Little h = H0 * 1e-2.
        mnu : float
            Sets the mass of a single neutrino species.
        E : float or array
            The normalised Hubble expansion rate.
        
        Returns
        -------
        Omega_nu : float or array
            Omega_nu fractional density.
        """
        return self.get_Omega_nu_E2(a, h, mnu)/(E**2)


    def compute_Omega_nu_prime(self, a, Omega_nu, E, E_prime, h, mnu):
        """
        Computes Omega photon prime.

        Parameters
        ----------
        Omega_g : float or array
            Radiation density.
        E : float or array
            Normalised Hubble expansion.
        E_prime : float or array
            Derivative of the normalised Hubble expansion.
        a : float or array
            Scale factor.
        h : float
            Little h = H0 * 1e-2.
        mnu : float
            Sets the mass of a single neutrino species.
        
        Returns
        -------
        Omega_nu_prime : float or array
            Omega neutrino prime.
        """
        E_prime_E = E_prime/E
        Omega_nu_prime = -(4 + 2*E_prime_E)*Omega_nu
        y = a * mnu / (self.params['Tnu0']*self.const['kB[eV]'])
        if np.isscalar(y):
            if y <= self.neutrino_table['ymin']:
                _Omega_nu_prime = self.neutrino_table['Iy_prime_interp_low'](y)
            elif y >= self.neutrino_table['ymax']:
                _Omega_nu_prime = self.neutrino_table['Iy_prime_interp_high'](y)
            else:
                _Omega_nu_prime = self.neutrino_table['Iy_prime_interp'](y)
        else:
            _Omega_nu_prime = np.zeros(len(y))
            cond = np.where((y <= self.neutrino_table['ymin']))[0]
            _Omega_nu_prime[cond] = self.neutrino_table['Iy_prime_interp_low'](y[cond])
            cond = np.where((y >= self.neutrino_table['ymax']))[0]
            _Omega_nu_prime[cond] = self.neutrino_table['Iy_prime_interp_high'](y[cond])
            cond = np.where((y > self.neutrino_table['ymin']) & (y < self.neutrino_table['ymax']))[0]
            _Omega_nu_prime[cond] = self.neutrino_table['Iy_prime_interp'](y[cond])
        _Omega_nu_prime *= self.const['C_nu']/((h**2) * (E**2) * (a**4))
        Omega_nu_prime += _Omega_nu_prime
        return Omega_nu_prime
    

    def compute_w_nu(self, a, mnu):
        """
        Computes the neutrino equation of state.

        Parameters
        ----------
        a : float or array
            Scale factor.
        mnu : float
            The mass of a single neutrino species.
        
        Returns
        -------
        w_nu : float or array
            Neutrino equation of state.
        """
        y = a * mnu / (self.params['Tnu0']*self.const['kB[eV]'])
        if np.isscalar(y):
            if y <= self.neutrino_table['ymin']:
                w_nu = self.neutrino_table['Jy_interp_low'](y)/(3*self.neutrino_table['Iy_interp_low'](y))
            elif y >= self.neutrino_table['ymax']:
                w_nu = self.neutrino_table['Jy_interp_high'](y)/(3*self.neutrino_table['Iy_interp_high'](y))
            else:
                w_nu = self.neutrino_table['Jy_interp'](y)/(3*self.neutrino_table['Iy_interp'](y))
        else:
            w_nu = np.zeros(len(y))
            cond = np.where((y <= self.neutrino_table['ymin']))[0]
            w_nu[cond] = self.neutrino_table['Jy_interp_low'](y[cond])/(3*self.neutrino_table['Iy_interp_low'](y[cond]))
            cond = np.where((y >= self.neutrino_table['ymax']))[0]
            w_nu[cond] = self.neutrino_table['Jy_interp_high'](y[cond])/(3*self.neutrino_table['Iy_interp_high'](y[cond]))
            cond = np.where((y > self.neutrino_table['ymin']) & (y < self.neutrino_table['ymax']))[0]
            w_nu[cond] = self.neutrino_table['Jy_interp'](y[cond])/(3*self.neutrino_table['Iy_interp'](y[cond]))
        return w_nu


    ## Matter
    
    def get_Omega_b_E2(self, a, Omega_b0):
        """
        Returns the baryon density.

        Parameters
        ----------
        a : float or array
            Scale factor.
        Omega_b0 : float
            Baryon fractional density at redshift zero.
        """
        return Omega_b0/(a**3)
    

    def get_Omega_b_E2_prime(self, a, Omega_b0):
        """
        Returns the baryon density derivative.

        Parameters
        ----------
        a : float or array
            Scale factor.
        Omega_b0 : float
            Baryon fractional density at redshift zero.
        """
        return -3*self.get_Omega_b_E2(a, Omega_b0)
    

    def get_Omega_b(self, a, Omega_b0, E):
        """
        Returns the baryon fractional density.

        Parameters
        ----------
        a : float or array
            Scale factor.
        Omega_b0 : float
            Baryon fractional density at redshift zero.
        E : float or array
            The normalised Hubble expansion rate.
        """
        return self.get_Omega_b_E2(a, Omega_b0)/(E**2)
    

    def compute_Omega_b_prime(self, Omega_b, E, E_prime):
        """
        Computes Omega baryon prime.

        Parameters
        ----------
        Omega_b : float or array
            Baryon density.
        E : float or array
            Normalised Hubble expansion.
        E_prime : float or array
            Derivative of the normalised Hubble expansion.
        
        Returns
        -------
        Omega_b_prime : float or array
            Omega baryon prime.
        """
        E_prime_E = E_prime/E
        Omega_b_prime = -Omega_b*(3. + 2.*E_prime_E)
        return Omega_b_prime


    def get_Omega_c_E2(self, a, Omega_c0):
        """
        Returns the cold dark matter density.

        Parameters
        ----------
        a : float or array
            Scale factor.
        Omega_c0 : float
            Cold dark matter fractional density at redshift zero.
        """
        return Omega_c0/(a**3)
    

    def get_Omega_c_E2_prime(self, a, Omega_c0):
        """
        Returns the cold dark matter density derivative.

        Parameters
        ----------
        a : float or array
            Scale factor.
        Omega_c0 : float
            Cold dark matter fractional density at redshift zero.
        """
        return -3*self.get_Omega_c_E2(a, Omega_c0)
    

    def get_Omega_c(self, a, Omega_c0, E):
        """
        Returns the cold dark matter fractional density.

        Parameters
        ----------
        a : float or array
            Scale factor.
        Omega_c0 : float
            Cold dark matter fractional density at redshift zero.
        E : float or array
            The normalised Hubble expansion rate.
        """
        return self.get_Omega_c_E2(a, Omega_c0)/(E**2)
    

    def compute_Omega_c_prime(self, Omega_c, E, E_prime):
        """
        Computes Omega cold dark matter prime.

        Parameters
        ----------
        Omega_b : float or array
            Cold dark matter density.
        E : float or array
            Normalised Hubble expansion.
        E_prime : float or array
            Derivative of the normalised Hubble expansion.
        
        Returns
        -------
        Omega_c_prime : float or array
            Omega cold dark matter prime.
        """
        E_prime_E = E_prime/E
        Omega_c_prime = -Omega_c*(3. + 2.*E_prime_E)
        return Omega_c_prime
    

    ## Dark energy

    def get_Omega_l_E2(self, a, Omega_l0, w0=-1., wa=0.):
        """
        Returns the dark energy density.

        Parameters
        ----------
        a : float or array
            Scale factor.
        Omega_l0 : float
            Baryon fractional density at redshift zero.
        w0 : float, optional
            Set dark energy equation of state today.
        wa : float, optional
            Set dark energy equation of state gradient.
        """
        if w0 == -1. and wa == 0.:
            return Omega_l0
        else:
            return Omega_l0*(a**(-3*(1+w0+wa)))*np.exp(3*wa*(a-1))
    

    def get_Omega_l_E2_prime(self, a, Omega_l0, w0=-1., wa=0.):
        """
        Returns the dark energy density.

        Parameters
        ----------
        a : float or array
            Scale factor.
        Omega_l0 : float
            Baryon fractional density at redshift zero.
        w0 : float, optional
            Set dark energy equation of state today.
        wa : float, optional
            Set dark energy equation of state gradient.
        """
        if w0 == -1. and wa == 0.:
            return 0.
        else:
            return -3*(1+w0+wa*(1-a))*self.get_Omega_l_E2(a, Omega_l0, w0=w0, wa=wa)
    

    def get_Omega_l(self, a, Omega_l0, E, w0=-1., wa=0.):
        """
        Returns the dark energy fractional density.

        Parameters
        ----------
        a : float or array
            Scale factor.
        Omega_l0 : float
            Baryon fractional density at redshift zero.
        E : float or array
            The normalised Hubble expansion rate.
        w0 : float, optional
            Set dark energy equation of state today.
        wa : float, optional
            Set dark energy equation of state gradient.
        """
        return self.get_Omega_l_E2(a, Omega_l0, w0=w0, wa=wa)/(E**2)
    

    def compute_Omega_l_prime(self, a, Omega_l, E, E_prime):
        """
        Computes Omega lambda prime.

        Parameters
        ----------
        Omega_l : float or array
            Cosmological constant energy density.
        E : float or array
            Normalised Hubble expansion.
        E_prime : float or array
            Derivative of the normalised Hubble expansion.
        
        Returns
        -------
        Omega_l_prime : float or array
            Omega Lambda or dynamical dark energy. 
        """
        E_prime_E = E_prime/E
        if self.params['w0'] == -1. and self.params['wa'] == 0.:
            Omega_l_prime = -2.*E_prime_E*Omega_l
        else:
            Omega_l_prime = -(2*E_prime_E + 3*(1+self.params['w0']+self.params['wa']*(1-a)))*Omega_l
            # Omega_l_prime = (-3*(1+self.params['w0']+self.params['wa']) + 3*self.params['wa']*a - 2.*E_prime_E)*Omega_l
        return Omega_l_prime
    

    def compute_w_l(self, a, w0=-1., wa=0.):
        """
        Returns the dark energy equation of state

        Parameters
        ----------
        a : float or array
            Scale factor.
        w0 : float, optional
            Set dark energy equation of state today.
        wa : float, optional
            Set dark energy equation of state gradient.
        """
        return w0 + wa*(1-a)

    
    # LCDM functions

    def compute_E_LCDM(self, a, h, Omega_b0, Omega_c0, mnu=[0.,0.,0.], w0=-1., wa=0.):
        """
        Computes the dimensionless Hubble expansion rate in LCDM.

        Parameters
        ----------
        a : float or array
            Scale factor.
        h : float
            Little h = H0 * 1e-2.
        Omega_b0 : float
            Baryon fractional density at redshift zero.
        Omega_c0 : float
            Cold dark matter fractional density at redshift zero.
        mnu : list, optional
            Sets the mass of the neutrino species.
        w0 : float, optional
            Set dark energy equation of state today.
        wa : float, optional
            Set dark energy equation of state gradient.

        Returns
        -------
        E : float or array
            The normalised Hubble expansion rate.
        """

        # Obtain photon and neutrino fractional density today to obtain dark energy density
        # assuming a spatially flat universe.

        Omega_g0 = self.get_Omega_g(1., h, 1.)
        Omega_nu0 = self.get_Omega_nu(1., h, mnu[0], 1.) + self.get_Omega_nu(1., h, mnu[1], 1.) + self.get_Omega_nu(1., h, mnu[2], 1.)
        
        Omega_l0 = 1. - Omega_g0 - Omega_nu0 - Omega_c0 - Omega_b0

        # Add each component to E2.

        E2 = self.get_Omega_g_E2(a, h)
        E2 += self.get_Omega_c_E2(a, Omega_c0)
        E2 += self.get_Omega_b_E2(a, Omega_b0)
        E2 += self.get_Omega_l_E2(a, Omega_l0, w0=w0, wa=wa)
        E2 += self.get_Omega_nu_E2(a, h, mnu[0]) + self.get_Omega_nu_E2(a, h, mnu[1]) + self.get_Omega_nu_E2(a, h, mnu[2])

        # Sqrt to obtain E

        E = np.sqrt(E2)

        return E
    

    def compute_E_prime_LCDM(self, a, h, Omega_c0, Omega_b0, mnu=[0.,0.,0.], w0=-1., wa=0.):
        """
        Computes the dimensionless Hubble expansion rate in LCDM.

        Parameters
        ----------
        a : float or array
            Scale factor.
        h : float
            Little h = H0 * 1e-2.
        Omega_c0 : float
            Cold dark matter fractional density at redshift zero.
        Omega_b0 : float
            Baryon fractional density at redshift zero.
        mnu : list, optional
            Sets the mass of the neutrino species.
        w0 : float, optional
            Set dark energy equation of state today.
        wa : float, optional
            Set dark energy equation of state gradient.

        Returns
        -------
        E_prime : float or array
            The derivative of normalised Hubble expansion rate.
        """

        E = self.compute_E_LCDM(a, h, Omega_c0, Omega_b0, mnu=mnu, w0=-1., wa=0.)

        # Obtain photon and neutrino fractional density today to obtain dark energy density
        # assuming a spatially flat universe.

        Omega_g0 = self.get_Omega_g(1., h, 1.)
        Omega_nu0 = self.get_Omega_nu(1., h, mnu[0], 1.) + self.get_Omega_nu(1., h, mnu[1], 1.) + self.get_Omega_nu(1., h, mnu[2], 1.)
        
        Omega_l0 = 1. - Omega_g0 - Omega_nu0 - Omega_c0 - Omega_b0

        E_prime = self.get_Omega_g_E2_prime(a, h)
        E_prime += self.get_Omega_nu_E2_prime(a, h, mnu[0]) + self.get_Omega_nu_E2_prime(a, h, mnu[1]) + self.get_Omega_nu_E2_prime(a, h, mnu[2])
        E_prime += self.get_Omega_b_E2_prime(a, Omega_b0)
        E_prime += self.get_Omega_c_E2_prime(a, Omega_c0)
        E_prime += self.get_Omega_l_E2_prime(a, Omega_l0, w0=w0, wa=wa)

        E_prime /= 2*E

        return E_prime


    def set_cosmo_params(self, H0_ref, Omega_c0_ref, Omega_b0_ref, fphi, K_G3_G4_values, w0=-1., wa=0., Tcmb=2.7255, Tnu0=1.9518, mnu=[0.,0.,0.]):
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
            Neutrino mass for each species.
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
        self.params['Omega_nu10_ref'] = self.get_Omega_nu(1., 1e-2*self.params['H0_ref'], self.params['mnu'][0], 1.) 
        self.params['Omega_nu20_ref'] = self.get_Omega_nu(1., 1e-2*self.params['H0_ref'], self.params['mnu'][1], 1.) 
        self.params['Omega_nu30_ref'] = self.get_Omega_nu(1., 1e-2*self.params['H0_ref'], self.params['mnu'][2], 1.)
        self.params['fphi'] = fphi
        self.params['Omega_l0_LCDM'] = 1. - self.params['Omega_g0_ref'] - self.params['Omega_nu10_ref'] - self.params['Omega_nu20_ref'] - self.params['Omega_nu30_ref']- self.params['Omega_c0_ref'] - self.params['Omega_b0_ref']
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

    # Returns the closure relation for a specific variable 

    def _fried_closure_wrapper(
            self, cl_val, cl_var, fried_closure_lambda, E_ini, phi_ini, phi_prime_ini, Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini, 
            Omega_nu1_ini, w_nu1_ini, Omega_nu2_ini, w_nu2_ini, Omega_nu3_ini, w_nu3_ini, w_l_ini, K_G3_G4_values, f_H_value):
        """
        Wrapper function for the Friedmann closure relation.

        Parameters
        ----------
        cl_val : float
            Closure variable value.
        cl_var : int
            Closure variable index, 0 = E, 1 = phi, 2 = phi_prime, 3 = Omega_g, 4 = Omega_b, 5 = Omega_c, 6 = Omega_l,
            7 = Omega_nu1, 8 = Omega_nu2, 9 = Omega_nu3.
        fried_closure_lamba : func
            Function for the closure relation, this should just `self.lambda_funcs['fried_closure_lambda']`.
        E_ini : float
            E initial value, if cl_var = 0 this value is computed through the closure relation.
        phi_ini : float
            phi initial value, if cl_var = 1 this value is computed through the closure relation.
        phi_prime_ini : float
            phi_prime initial value, if cl_var = 2 this value is computed through the closure relation.
        Omega_g_ini : float
            Omega photon initial value, if cl_var = 3 this value is computed through the closure relation.
        Omega_b_prime_ini : float
            Omega baryon initial value, if cl_var = 4 this value is computed through the closure relation.
        Omega_c_prime_ini : float
            Omega cold dark matter initial value, if cl_var = 5 this value is computed through the closure relation.
        Omega_l_ini : float
            Omega lambda initial value, if cl_var = 6 this value is computed through the closure relation.
        Omega_nu1_ini : float
            Omega neutrino species 1 initial value, if cl_var = 7 this value is computed through the closure relation.
        w_nu1_ini : float
            Neutrino species 1 equation of state.
        Omega_nu2_ini : float
            Omega neutrino species 2 initial value, if cl_var = 8 this value is computed through the closure relation.
        w_nu2_ini : float
            Neutrino species 2 equation of state.
        Omega_nu3_ini : float
            Omega neutrino species 3 initial value, if cl_var = 9 this value is computed through the closure relation.
        w_nu3_ini : float
            Neutrino species 3 equation of state.
        w_l_ini : float
            Dark energy equation of state.
        K_G3_G4_values : float
            A list of values for the Horndeski specific variables. This must match the length of the user defined variable. 
            Check self.sym['K_G3_G4_syms'] to see what variables are expected.
        """
        assert cl_var >= 0 and cl_var <= 9, "Closure variable unsupported, must be between 0 and 9 inclusive."
        if cl_var == 0:
            #Closure used to compute E0
            return fried_closure_lambda(cl_val, phi_ini, phi_prime_ini, Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini, 
                                        Omega_nu1_ini, w_nu1_ini, Omega_nu2_ini, w_nu2_ini, Omega_nu3_ini, w_nu3_ini, w_l_ini, *K_G3_G4_values, f_H_value)
        if cl_var == 1:
            #Closure used to compute phi0
            return fried_closure_lambda(E_ini, cl_val, phi_prime_ini, Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini, 
                                        Omega_nu1_ini, w_nu1_ini, Omega_nu2_ini, w_nu2_ini, Omega_nu3_ini, w_nu3_ini, w_l_ini, *K_G3_G4_values, f_H_value)
        if cl_var == 2:
            #Closure used to compute phi_prime0
            return fried_closure_lambda(E_ini, phi_ini, cl_val, Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini, 
                                        Omega_nu1_ini, w_nu1_ini, Omega_nu2_ini, w_nu2_ini, Omega_nu3_ini, w_nu3_ini, w_l_ini, *K_G3_G4_values, f_H_value)
        if cl_var == 3:
            #Closure used to compute Omega_g0
            return fried_closure_lambda(E_ini, phi_ini, phi_prime_ini, cl_val, Omega_b_ini, Omega_c_ini, Omega_l_ini, 
                                        Omega_nu1_ini, w_nu1_ini, Omega_nu2_ini, w_nu2_ini, Omega_nu3_ini, w_nu3_ini, w_l_ini, *K_G3_G4_values, f_H_value)
        if cl_var == 4:
            #Closure used to compute Omega_b0
            return fried_closure_lambda(E_ini, phi_ini, phi_prime_ini, Omega_g_ini, cl_val, Omega_c_ini, Omega_l_ini, 
                                        Omega_nu1_ini, w_nu1_ini, Omega_nu2_ini, w_nu2_ini, Omega_nu3_ini, w_nu3_ini, w_l_ini, *K_G3_G4_values, f_H_value)
        if cl_var == 5:
            #Closure used to compute Omega_c0
            return fried_closure_lambda(E_ini, phi_ini, phi_prime_ini, Omega_g_ini, Omega_b_ini, cl_val, Omega_l_ini, 
                                        Omega_nu1_ini, w_nu1_ini, Omega_nu2_ini, w_nu2_ini, Omega_nu3_ini, w_nu3_ini, w_l_ini, *K_G3_G4_values, f_H_value)
        if cl_var == 6:
            #Closure used to compute Omega_l0
            return fried_closure_lambda(E_ini, phi_ini, phi_prime_ini, Omega_g_ini, Omega_b_ini, Omega_c_ini, cl_val, 
                                        Omega_nu1_ini, w_nu1_ini, Omega_nu2_ini, w_nu2_ini, Omega_nu3_ini, w_nu3_ini, w_l_ini, *K_G3_G4_values, f_H_value)
        if cl_var == 7:
            #Closure used to compute Omega_nu10
            return fried_closure_lambda(E_ini, phi_ini, phi_prime_ini, Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini, 
                                        Omega_nu1_ini, w_nu1_ini, Omega_nu2_ini, w_nu2_ini, Omega_nu3_ini, w_nu3_ini, w_l_ini, *K_G3_G4_values, f_H_value)
        if cl_var == 8:
            #Closure used to compute Omega_nu20
            return fried_closure_lambda(E_ini, phi_ini, phi_prime_ini, Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini, 
                                        cl_val, w_nu1_ini, Omega_nu2_ini, w_nu2_ini, Omega_nu3_ini, w_nu3_ini, w_l_ini, *K_G3_G4_values, f_H_value)
        if cl_var == 9:
            #Closure used to compute Omega_nu30
            return fried_closure_lambda(E_ini, phi_ini, phi_prime_ini, Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini, 
                                        Omega_nu1_ini, w_nu1_ini, cl_val, w_nu2_ini, cl_val, w_nu3_ini, w_l_ini, *K_G3_G4_values, f_H_value)


    # Computes derivaties for the solver

    def _compute_primes(self, x, Y, K_G3_G4_values, f_H_value, timeout=5):
        """
        Compute prime functions for numerical solver.

        Parameters
        ----------
        x : float
            Current value of log(a).
        Y : list
            List containing current [_E, _phi, _phi_prime, _Omega_g, _Omega_b, _Omega_c, _Omega_l, _Omega_nu1, _Omega_nu2, _Omega_nu3] values.
        K_G3_G4_values : list
            A list of values for the Horndeski specific variables. This must match the length of the user defined variable. 
            Check self.sym['K_G3_G4_syms'] to see what variables are expected.
        timeout : float, optional
            Time in seconds to force the solver to fail.
        """
        # convert x = log(a) to scale factor
        a = np.exp(x)

        # `_` used to denote current value.
        _E, _phi, _phi_prime, _Omega_g, _Omega_b, _Omega_c, _Omega_l, _Omega_nu1, _Omega_nu2, _Omega_nu3 = Y

        _w_nu1 = self.compute_w_nu(a, self.params['mnu'][0])
        _w_nu2 = self.compute_w_nu(a, self.params['mnu'][1])
        _w_nu3 = self.compute_w_nu(a, self.params['mnu'][2])
        _w_l = self.compute_w_l(a, w0=self.params['w0'], wa=self.params['wa'])

        variables = [_E, _phi, _phi_prime, _Omega_g, _Omega_b, _Omega_c, _Omega_l, 
                     _Omega_nu1, _w_nu1, _Omega_nu2, _w_nu2, _Omega_nu3, _w_nu3, _w_l, *K_G3_G4_values, f_H_value]

        phi_prime = _phi_prime
        E_prime_E = self.lambda_funcs['E_prime/E_lambda'](*variables)
        E_prime = E_prime_E*_E
        phi_primeprime = self.lambda_funcs['phi_primeprime_lambda'](E_prime, *variables)

        Omega_g_prime = self.compute_Omega_g_prime(_Omega_g, _E, E_prime)
        Omega_b_prime = self.compute_Omega_b_prime(_Omega_b, _E, E_prime)
        Omega_c_prime = self.compute_Omega_c_prime(_Omega_c, _E, E_prime)
        Omega_l_prime = self.compute_Omega_l_prime(a, _Omega_l, _E, E_prime)
        Omega_nu1_prime = self.compute_Omega_nu_prime(a, _Omega_nu1, _E, E_prime, 1e-2*self.params['H0_ref'], self.params['mnu'][0])
        Omega_nu2_prime = self.compute_Omega_nu_prime(a, _Omega_nu2, _E, E_prime, 1e-2*self.params['H0_ref'], self.params['mnu'][1])
        Omega_nu3_prime = self.compute_Omega_nu_prime(a, _Omega_nu3, _E, E_prime, 1e-2*self.params['H0_ref'], self.params['mnu'][2])

        timenow = self._check_timer()

        if timenow >= timeout:
            
            E_prime = np.nan
            phi_prime = np.nan
            phi_primeprime = np.nan
            Omega_g_prime = np.nan
            Omega_b_prime = np.nan
            Omega_c_prime = np.nan
            Omega_l_prime = np.nan
            Omega_nu1_prime = np.nan
            Omega_nu2_prime = np.nan
            Omega_nu3_prime = np.nan

        if np.isfinite([
            E_prime, phi_prime, phi_primeprime, Omega_g_prime, 
            Omega_b_prime, Omega_c_prime, Omega_l_prime, 
            Omega_nu1_prime, Omega_nu2_prime, Omega_nu3_prime]).all() == False:
            self._solver_success = False

        return [E_prime, phi_prime, phi_primeprime, Omega_g_prime, Omega_b_prime, Omega_c_prime, Omega_l_prime, Omega_nu1_prime, Omega_nu2_prime, Omega_nu3_prime]

    # Numerically computed quantities

    def compute_chi_over_delta(self, a, E, calB, calC):
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

        Return
        ------
        chioverdelta : float or array   
            Code equivalent of equation 3.14 in https://arxiv.org/abs/2209.01666.
        """
        if self.params['Omega_c0'] is None:
            chioverdelta = None
        else:
            chioverdelta = calB * calC * (self.params['Omega_c0']+self.params['Omega_b0'])/((E**2)*(a**3)) # TODO: there's a G_G4/G_N in 3.14 in https://arxiv.org/pdf/2209.01666 which is not included here...
        return chioverdelta
    

    def comp_w_eff(self):
        """
        Computes effective equation of state.
        """
        if self.output['E_prime/E'] is None:
            self.output['w_eff'] = None
        else:
            self.output['w_eff'] = - 1 - (2/3)*self.output['E_prime/E']
    

    def get_phi_sound_horizon(self):
        """
        Computes the sound horizon for the scalar field.
        """

        if self.output['E'] is None:

            r_s = None
            r_s_c1 = None

        else:
            
            from scipy.integrate import cumulative_trapezoid

            if self.output['E'].ndim == 1:
                
                keys = ['E', 'c_s_sq']

                check = True
                for key in keys:
                    if self.output[key] is None or np.isfinite(self.output[key]).all() == False:
                        check = False

                if check:
                    c_s = np.sqrt(self.output['c_s_sq'])*self.const['c[km/s]']
                    integral = c_s / (100.*np.exp(self.output['x']) * self.output['E'])
                    r_s = cumulative_trapezoid(integral, x=self.output['x'], initial=0.)
                    c_s1 = np.ones(len(self.output['c_s_sq']))*self.const['c[km/s]']
                    integral = c_s1 / (100.*np.exp(self.output['x']) * self.output['E'])
                    r_s_c1 = cumulative_trapezoid(integral, x=self.output['x'], initial=0.)
                else:
                    r_s = None
                    r_s_c1 = None
        
            elif self.output['E'].ndim == 2:
        
                r_s = np.zeros(np.shape(self.output['E']))
                r_s_c1 = np.zeros(np.shape(self.output['E']))

                for idx in range(0, len(self.output['Omega_c0'])):
                    
                    keys = ['E', 'c_s_sq']

                    check = True
                    for key in keys:
                        if self.output[key] is None or np.isfinite(self.output[key][idx]).all() == False:
                            check = False

                    if check:
                        c_s = np.sqrt(self.output['c_s_sq'][idx])*self.const['c[km/s]']
                        integral = c_s / (100.*np.exp(self.output['x']) * self.output['E'][idx])
                        r_s[idx] = cumulative_trapezoid(integral, x=self.output['x'], initial=0.)
                        c_s1 = np.ones(len(self.output['c_s_sq'][idx]))*self.const['c[km/s]']
                        integral = c_s1 / (100.*np.exp(self.output['x']) * self.output['E'][idx])
                        r_s_c1[idx] = cumulative_trapezoid(integral, x=self.output['x'], initial=0.)
                    else:
                        r_s[idx] = np.nan * np.ones(len(self.output['x']))
                        r_s_c1[idx] = np.nan * np.ones(len(self.output['x']))
            else:

                r_s = None
                r_s_c1 = None

        self.output['r_s'] = r_s
        self.output['r_s_c1'] = r_s_c1


    def _linear_growth(self, x, y, Omega_m0):
        """
        Linear growth ODE system.

        Parameters
        ----------
        x : float
            Log of the scale factor.
        y : float
            The growth function and it's derivative [D, dD].
        Omega_m0 : float
            The matter density.
        
        Returns
        -------
        dD : float
            Derivative of the growth function.
        d2D : float
            The second order derivative of the growth function.
        """
        D = y[0]
        dD = y[1]
        a = np.exp(x)
        mu = 1 + self._linear_growth_int['interp_beta_vs_a'](a)
        GG4GN = self._linear_growth_int['interp_G_G_4/G_N_vs_a'](a)
        Bx = 2 + self._linear_growth_int['interp_E_prime_vs_a'](a)/self._linear_growth_int['interp_E_vs_a'](a)
        Cx = 3.*Omega_m0*mu*GG4GN/(2.*(self._linear_growth_int['interp_E_vs_a'](a)**2)*(a**3))
        d2D = Cx*D - Bx*D
        return [dD, d2D]
    

    def get_linear_growth(self):
        """
        Compute the first order linear growth function.
        """
        
        from scipy.interpolate import interp1d
        from scipy.integrate import solve_ivp

        # Run sanity checks to test whether linear growth functions can be computed
        if self.output['E'] is None:
            
            D1, f1 = None, None

        else:

            if self.output['E'].ndim == 1:
                
                keys = ['E', 'G_G_4/G_N', 'E_prime', 'beta']

                check = True

                for key in keys:
                    if self.output[key] is None or np.isfinite(self.output[key]).all() == False:
                        check = False
                
                if check:

                    self._linear_growth_int = {}
                    self._linear_growth_int['interp_E_vs_a'] = interp1d(self.output['a'], self.output['E'], kind='cubic', fill_value='extrapolate')
                    self._linear_growth_int['interp_E_prime_vs_a'] = interp1d(self.output['a'], self.output['E_prime'], kind='cubic', fill_value='extrapolate')
                    self._linear_growth_int['interp_beta_vs_a'] = interp1d(self.output['a'], self.output['beta'], kind='cubic', fill_value='extrapolate')
                    self._linear_growth_int['interp_G_G_4/G_N_vs_a'] = interp1d(self.output['a'], self.output['G_G_4/G_N'], kind='cubic', fill_value='extrapolate')

                    # set up initial conditons, assuming matter domination.
                    
                    x_ini = self.output['x'][0]
                    D_ini = self.output['a'][0]
                    dD_ini = self.output['a'][0]
                    y_ini = [D_ini, dD_ini]

                    # End position
                    x_final = self.output['x'][-1]
                    
                    # Solve forward
                    ans = solve_ivp(self._linear_growth, (x_ini, x_final), y_ini, t_eval=self.output['x'], args=(self.output['Omega_c0']+self.output['Omega_b0'],))
                        
                    # Combine solutions
                    _D1 = ans.y[0]
                    _dD1 = ans.y[1]

                    _D1_constant = np.copy(_D1[-1])
                    _D1 /= _D1_constant
                    _D1prime = _dD1
                    _D1prime /= _D1_constant

                    _f1 = _D1prime/_D1

                    D1, f1 = _D1, _f1
                
                else:

                    D1, f1 = None, None
                
            elif self.output['E'].ndim == 2:
            
                D1 = np.zeros(np.shape(self.output['E']))
                f1 = np.zeros(np.shape(self.output['E']))

                for idx in range(0, len(self.output['Omega_c0'])):
                    
                    keys = ['E', 'G_G_4/G_N', 'E_prime', 'beta']

                    check = True

                    for key in keys:
                        if self.output[key] is None or np.isfinite(self.output[key][idx]).all() == False:
                            check = False
                    
                    if check:

                        self._linear_growth_int = {}
                        self._linear_growth_int['interp_E_vs_a'] = interp1d(self.output['a'], self.output['E'][idx], kind='cubic', fill_value='extrapolate')
                        self._linear_growth_int['interp_E_prime_vs_a'] = interp1d(self.output['a'], self.output['E_prime'][idx], kind='cubic', fill_value='extrapolate')
                        self._linear_growth_int['interp_beta_vs_a'] = interp1d(self.output['a'], self.output['beta'][idx], kind='cubic', fill_value='extrapolate')
                        self._linear_growth_int['interp_G_G_4/G_N_vs_a'] = interp1d(self.output['a'], self.output['G_G_4/G_N'][idx], kind='cubic', fill_value='extrapolate')

                        # set up initial conditons, assuming matter domination.
                        
                        x_ini = self.output['x'][0]
                        D_ini = self.output['a'][0]
                        dD_ini = self.output['a'][0]
                        y_ini = [D_ini, dD_ini]

                        # End position
                        x_final = self.output['x'][-1]
                        
                        # Solve forward
                        ans = solve_ivp(self._linear_growth, (x_ini, x_final), y_ini, t_eval=self.output['x'], args=(self.output['Omega_c0'][idx]+self.output['Omega_b0'][idx],))
                            
                        # Combine solutions
                        _D1 = ans.y[0]
                        _dD1 = ans.y[1]

                        _D1_constant = np.copy(_D1[-1])
                        _D1 /= _D1_constant
                        _D1prime = _dD1
                        _D1prime /= _D1_constant

                        _f1 = _D1prime/_D1

                        D1[idx], f1[idx] = _D1, _f1

                    else:
                        D1[idx] = np.nan * np.ones(len(self.output['x']))
                        f1[idx] = np.nan * np.ones(len(self.output['x']))

            else:
                
                D1, f1 = None, None

        self.output['D1'] = D1
        self.output['f1'] = f1

    
    def _linear_growth_2(self, x, y, Omega_m0):
        """
        Linear growth ODE system.

        Parameters
        ----------
        x : float
            Log of the scale factor.
        y : float
            The growth function and it's derivative [D, dD].
        Omega_m0 : float
            The matter density.
        
        Returns
        -------
        dD2 : float
            Derivative of the second order growth function.
        d2D2 : float
            The second order derivative of the second order growth function.
        """
        D2 = y[0]
        dD2 = y[1]
        a = np.exp(x)
        D1 = self._linear_growth_int['interp_D1_vs_a'](a)
        # assuming mu2 == mu1
        mu = 1 + self._linear_growth_int['interp_beta_vs_a'](a)
        GG4GN = self._linear_growth_int['interp_G_G_4/G_N_vs_a'](a)
        Bx = 2 + self._linear_growth_int['interp_E_prime_vs_a'](a)/self._linear_growth_int['interp_E_vs_a'](a)
        Cx = 3.*Omega_m0*mu*GG4GN/(2.*(self._linear_growth_int['interp_E_vs_a'](a)**2)*(a**3))
        d2D2 = Cx*(D2 - D1**2) - Bx*dD2
        return [dD2, d2D2]


    def get_linear_growth_2(self):
        """
        Compute the second order linear growth function assuming mu1 = mu2.
        """
        
        from scipy.interpolate import interp1d
        from scipy.integrate import solve_ivp

        # Run sanity checks to test whether linear growth functions can be computed
        if self.output['E'] is None:
            
            D2, f2 = None, None
        
        else:

            if self.output['E'].ndim == 1:
                
                keys = ['E', 'E_prime', 'beta', 'D1']
                
                check = True
                for key in keys:

                    if self.output[key] is None or np.isfinite(self.output[key]).all() == False:
                        check = False

                if check:
                        
                    self._linear_growth_int = {}
                    self._linear_growth_int['interp_E_vs_a'] = interp1d(self.output['a'], self.output['E'], kind='cubic', fill_value='extrapolate')
                    self._linear_growth_int['interp_E_prime_vs_a'] = interp1d(self.output['a'], self.output['E_prime'], kind='cubic', fill_value='extrapolate')
                    self._linear_growth_int['interp_G_G_4/G_N_vs_a'] = interp1d(self.output['a'], self.output['G_G_4/G_N'], kind='cubic', fill_value='extrapolate')
                    self._linear_growth_int['interp_beta_vs_a'] = interp1d(self.output['a'], self.output['beta'], kind='cubic', fill_value='extrapolate')
                    self._linear_growth_int['interp_D1_vs_a'] = interp1d(self.output['a'], self.output['D1'], kind='cubic', fill_value='extrapolate')

                    # set up initial conditons, assuming matter domination.
                    
                    x_ini = self.output['x'][0]
                    D2_ini = -(3/7)*self.output['a'][0]**2
                    dD2_ini = -(6/7)*self.output['a'][0]**2
                    
                    y_ini = [D2_ini, dD2_ini]

                    # End position
                    x_final = self.output['x'][-1]
                    
                    # Solve forward
                    ans = solve_ivp(self._linear_growth_2, (x_ini, x_final), y_ini, t_eval=self.output['x'], args=(self.output['Omega_c0']+self.output['Omega_b0'],))
                        
                    # Combine solutions
                    _D2 = ans.y[0]
                    _dD2 = ans.y[1]

                    _f2 = _dD2/_D2

                    D2, f2 = _D2, _f2
            
                else:

                    D2, f2 = None, None

            elif self.output['E'].ndim == 2:
                
                D2 = np.zeros(np.shape(self.output['E']))
                f2 = np.zeros(np.shape(self.output['E']))

                for idx in range(0, len(self.output['Omega_c0'])):
                    
                    keys = ['E', 'E_prime', 'beta', 'D1']
                    
                    check = True
                    for key in keys:

                        if self.output[key] is None or np.isfinite(self.output[key][idx]).all() == False:
                            check = False

                    if check:

                        self._linear_growth_int = {}
                        self._linear_growth_int['interp_E_vs_a'] = interp1d(self.output['a'], self.output['E'][idx], kind='cubic', fill_value='extrapolate')
                        self._linear_growth_int['interp_E_prime_vs_a'] = interp1d(self.output['a'], self.output['E_prime'][idx], kind='cubic', fill_value='extrapolate')
                        self._linear_growth_int['interp_G_G_4/G_N_vs_a'] = interp1d(self.output['a'], self.output['G_G_4/G_N'][idx], kind='cubic', fill_value='extrapolate')
                        self._linear_growth_int['interp_beta_vs_a'] = interp1d(self.output['a'], self.output['beta'][idx], kind='cubic', fill_value='extrapolate')
                        self._linear_growth_int['interp_D1_vs_a'] = interp1d(self.output['a'], self.output['D1'][idx], kind='cubic', fill_value='extrapolate')

                        # set up initial conditons, assuming matter domination.
                        
                        x_ini = self.output['x'][0]
                        D2_ini = -(3/7)*self.output['a'][0]**2
                        dD2_ini = -(6/7)*self.output['a'][0]**2

                        y_ini = [D2_ini, dD2_ini]

                        # End position
                        x_final = self.output['x'][-1]
                        
                        # Solve forward
                        ans = solve_ivp(self._linear_growth_2, (x_ini, x_final), y_ini, t_eval=self.output['x'], args=(self.output['Omega_c0'][idx]+self.output['Omega_b0'][idx],))
                            
                        # Combine solutions
                        _D2 = ans.y[0]
                        _dD2 = ans.y[1]

                        _f2 = _dD2/_D2

                        D2[idx], f2[idx] = _D2, _f2
                
                    else:
                        
                        D2[idx] = np.nan * np.ones(len(self.output['x']))
                        f2[idx] = np.nan * np.ones(len(self.output['x']))

            else:

                D2, f2 = None, None
        
        self.output['D2'] = D2
        self.output['f2'] = f2


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
                        beta1 = 3*(self.output['rho_phi']/(self.output['H']**2))*(1 + self.output['w_phi']) + self.output['E_prime/E']*(4 + self.output['alpha_B']) + self.output['alpha_B_prime']
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
                            beta1 = 3*(self.output['rho_phi'][idx]/(self.output['H'][idx]**2))*(1 + self.output['w_phi'][idx]) + self.output['E_prime/E'][idx]*(4 + self.output['alpha_B'][idx]) + self.output['alpha_B_prime'][idx]
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
    
        
    def get_Sigma_derivatives(self):
        """
        Computes the derivates of the Sigma useful for computing the ISW.
        """

        if self.output['E'] is None:

            Sigma1, S, zeta = None, None, None
        
        else:

            if self.output['E'].ndim == 1:

                check = True

                if self.output['Sigma'] is None or np.isfinite(self.output['Sigma']).all() == False:
                    check = False

                if check:
                    Sigma1 = self.output['Sigma'][-1]
                    S = self.output['Sigma']/Sigma1
                    zeta = np.gradient(np.log(S), self.output['x'])
                else:
                    Sigma1, S, zeta = None, None, None
                
            elif self.output['E'].ndim == 2:

                Sigma1 = np.zeros(np.shape(self.output['E']))
                S = np.zeros(np.shape(self.output['E']))
                zeta = np.zeros(np.shape(self.output['E']))

                for idx in range(0, len(self.output['Omega_c0'])):

                    check = True

                    if self.output['beta'][idx] is None or np.isfinite(self.output['beta'][idx]).all() == False:
                        check = False

                    if check:
                        Sigma1[idx] = self.output['Sigma'][idx][-1]
                        S[idx] = self.output['Sigma'][idx]/Sigma1[idx]
                        zeta[idx] = np.gradient(np.log(S[idx]), self.output['x'])
                    else:
                        Sigma1[idx] = np.nan
                        S[idx] = np.nan * np.ones(len(self.output['x']))
                        zeta[idx] = np.nan * np.ones(len(self.output['x']))
            else:

                Sigma1, S, zeta = None, None, None
        
        self.output['Sigma1'] = Sigma1
        self.output['S'] = S
        self.output['zeta'] = zeta


    # alternative numerical computation of w_phi

    def compute_w_phi(self):
        """
        Computes equation of state for scalar field using the Friedman equations. 
        Code equivalent of 3.5 in https://iopscience.iop.org/article/10.1088/1475-7516/2014/07/050.
        """
        # TODO : need to look at how this is computed.
        if self.output['E'] is None:

            P_phi, rho_phi, w_phi = None, None, None

        else:
            
            if self.output['E'].ndim == 1:
                
                check = True
                keys = ['E', 'E_prime', 'Omega_g', 'Omega_c', 'Omega_b', 'Omega_l', 'Omega_nu1', 'Omega_nu2', 'Omega_nu3', 'M_star_sq']

                for key in keys:
                    if np.isfinite(self.output[key]).all() == False or self.output[key] is None:
                        check = False 
                
                if check:

                    P_tot = self.output['Omega_g']/3.
                    P_tot += self.compute_w_l(self.output['a'], self.params['w0'], self.params['wa'])*self.output['Omega_l']
                    # P_tot += self.compute_w_nu(self.output['a'], self.params['mnu'][0])*self.output['Omega_nu1']
                    # P_tot += self.compute_w_nu(self.output['a'], self.params['mnu'][1])*self.output['Omega_nu2']
                    # P_tot += self.compute_w_nu(self.output['a'], self.params['mnu'][2])*self.output['Omega_nu3']

                    P_tot += self.output['w_nu1']*self.output['Omega_nu1']
                    P_tot += self.output['w_nu2']*self.output['Omega_nu2']
                    P_tot += self.output['w_nu3']*self.output['Omega_nu3']

                    rho_tot = self.output['Omega_g'] + self.output['Omega_b'] + self.output['Omega_c'] + self.output['Omega_l']
                    rho_tot += self.output['Omega_nu1'] + self.output['Omega_nu2'] + self.output['Omega_nu3']

                    P_phi = -2*(self.params['H0']**2)*self.output['E']*self.output['E_prime'] 
                    P_phi -= 3*(self.params['H0']**2)*(self.output['E']**2)*(1+P_tot/self.output['M_star_sq'])
                    rho_phi = 3*(self.params['H0']**2)*(self.output['E']**2)*(1-rho_tot/self.output['M_star_sq'])

                    w_phi = np.zeros(len(self.output['rho_phi']))
                    cond = np.where(self.output['rho_phi'] == 0.)[0]
                    w_phi[cond] = np.nan
                    cond = np.where(self.output['rho_phi'] != 0.)[0]
                    w_phi[cond] = self.output['P_phi'][cond]/self.output['rho_phi'][cond]
                else:
                    P_phi, rho_phi, w_phi = None, None, None

            elif self.output['E'].ndim == 2:
                
                P_phi = np.zeros(np.shape(self.output['E']))
                rho_phi = np.zeros(np.shape(self.output['E']))
                w_phi = np.zeros(np.shape(self.output['E']))

                for idx in range(0, len(self.output['Omega_c0'])):
                    
                    check = True
                    keys = ['E', 'E_prime', 'Omega_g', 'Omega_c', 'Omega_b', 'Omega_l', 'Omega_nu1', 'Omega_nu2', 'Omega_nu3', 'M_star_sq']

                    for key in keys:
                        if np.isfinite(self.output[key][idx]).all() == False or self.output[key] is None:
                            check = False 

                    if check:
                        P_tot = self.output['Omega_g'][idx]/3.
                        P_tot += self.compute_w_l(self.output['a'], self.params['w0'], self.params['wa'])*self.output['Omega_l'][idx]
                        P_tot += self.output['w_nu1'][idx]*self.output['Omega_nu1'][idx]
                        P_tot += self.output['w_nu2'][idx]*self.output['Omega_nu2'][idx]
                        P_tot += self.output['w_nu3'][idx]*self.output['Omega_nu3'][idx]

                        rho_tot = self.output['Omega_g'][idx] + self.output['Omega_b'][idx] + self.output['Omega_c'][idx] + self.output['Omega_l'][idx]
                        rho_tot += self.output['Omega_nu1'][idx] + self.output['Omega_nu2'][idx] + self.output['Omega_nu3'][idx]


                        P_phi[idx] = -2*(self.params['H0']**2)*self.output['E'][idx]*self.output['E_prime'][idx] 
                        P_phi[idx] -= 3*(self.params['H0']**2)*(self.output['E'][idx]**2)*(1+P_tot/self.output['M_star_sq'][idx])
                        rho_phi[idx] = 3*(self.params['H0']**2)*(self.output['E'][idx]**2)*(1-rho_tot/self.output['M_star_sq'][idx])
                        w_phi[idx] = P_phi[idx]/rho_phi[idx]
                    else:
                        P_phi[idx] = np.nan * np.ones(len(self.output['x']))
                        rho_phi[idx] = np.nan * np.ones(len(self.output['x']))
                        w_phi[idx] = np.nan * np.ones(len(self.output['x']))

            else:

                P_phi, rho_phi, w_phi = None, None, None
        
        self.output['P_phi_alt'] = P_phi
        self.output['rho_phi_alt'] = rho_phi
        self.output['w_phi_alt'] = w_phi
    

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

        E_ini = self.compute_E_LCDM(a_start, 1e-2*self.params['H0_ref'], self.params['Omega_b0_ref'], self.params['Omega_c0_ref'], mnu=self.params['mnu'], w0=self.params['w0'], wa=self.params['wa'])
        #lcdm.compute_Ez_LCDM(z_start, self.params['Omega_r0_ref'], self.params['Omega_m0_ref'], w0=self.params['w0'], wa=self.params['wa'])
        Omega_g_ini = self.get_Omega_g(a_start, 1e-2*self.params['H0_ref'], E_ini)
        #Omega_r_ini = lcdm.compute_Omega_r_z_LCDM(z_start, self.params['Omega_r0_ref'], self.params['Omega_m0_ref'], w0=self.params['w0'], wa=self.params['wa'])
        Omega_nu1_ini = self.get_Omega_nu(a_start, 1e-2*self.params['H0_ref'], self.params['mnu'][0], E_ini)
        Omega_nu2_ini = self.get_Omega_nu(a_start, 1e-2*self.params['H0_ref'], self.params['mnu'][1], E_ini)
        Omega_nu3_ini = self.get_Omega_nu(a_start, 1e-2*self.params['H0_ref'], self.params['mnu'][2], E_ini)
        
        Omega_b_ini = self.get_Omega_b(a_start, self.params['Omega_b0_ref'], E_ini)
        Omega_c_ini = self.get_Omega_c(a_start, self.params['Omega_c0_ref'], E_ini)
        #Omega_m_ini = lcdm.compute_Omega_m_z_LCDM(z_start, self.params['Omega_r0_ref'], self.params['Omega_m0_ref'], w0=self.params['w0'], wa=self.params['wa'])

        Omega_l_ini = (1.-self.params['fphi'])*self.get_Omega_l(a_start, self.params['Omega_l0_LCDM'], E_ini, w0=self.params['w0'], wa=self.params['wa'])
        #Omega_l_ini = (1.-self.params['fphi'])*lcdm.compute_Omega_l_z_LCDM(z_start, self.params['Omega_r0_ref'], self.params['Omega_m0_ref'], w0=self.params['w0'], wa=self.params['wa'])

        return E_ini, Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini, Omega_nu1_ini, Omega_nu2_ini, Omega_nu3_ini
    

    def _run_solver_HG_root_finder(
            self, closure_variable, E_ini, phi_ini, phi_prime_ini, Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini, 
            Omega_nu1_ini, Omega_nu2_ini, Omega_nu3_ini, bypass_closure):
        """
        Finds the roots of the closure relation for input closure variable given a set of initial guesses. Note the initial guess for the 
        variable of interest will be ignored.

        Parameters
        ----------
        closure_variable : str
            Variable used to set the initial conditions.
        E_ini : float
            Initial normalised Hubble factor, if the closure_variable = 0 then this is ignored and solved via the closure relation.
        phi_ini : float
            phi initial value, if the closure_variable = 1 then this is ignored and solved via the closure relation.
        phi_prime_ini : float
            phi_prime initial value, if the closure_variable = 2 then this is ignored and solved via the closure relation.
        Omega_g_ini : float
            Initial fractional radiation density, if the closure_variable = 3 then this is ignored and solved via the closure relation.
        Omega_b_ini : float
            Initial fractional baryon density, if the closure_variable = 4 then this is ignored and solved via the closure relation.
        Omega_c_ini : float
            Initial fractional cold dark matter density, if the closure_variable = 5 then this is ignored and solved via the closure relation.
        Omega_l_ini : float
            Initial fractional lambda (cosmological constant) density, if the closure_variable = 6 then this is ignored and solved via the closure relation.
        Omega_nu1_ini : float
            Initial neutrino density for neutrino species 1, if closure_variable = 7 then this is ignored.
        Omega_nu2_ini : float
            Initial neutrino density for neutrino species 2, if closure_variable = 8 then this is ignored.
        Omega_nu3_ini : float
            Initial neutrino density for neutrino species 3, if closure_variable = 9 then this is ignored.
        bypass_closure : bool, optional
            Bypass closure solver to set quantities directly, this means the closure relation will not be satisfied, use with care!
        """
        
        sub_dict = {}
            
        assert closure_variable >= 0 and closure_variable <= 5, "Closure variable unsupported, must be between 0 and 5 inclusive."

        if closure_variable != 0:
            sub_dict[self.sym['E']] = E_ini
        else:
            closure_string = 'E'
            closure_variable_sym = self.sym['E']
        if closure_variable != 1:
            sub_dict[self.sym['phi']] = phi_ini
        else:
            closure_string = 'phi'
            closure_variable_sym = self.sym['phi']
        if closure_variable != 2:
            sub_dict[self.sym['phi_prime']] = phi_prime_ini
        else:
            closure_string = 'phi_prime'
            closure_variable_sym = self.sym['phi_prime']
        if closure_variable != 3:
            sub_dict[self.sym['Omega_g']] = Omega_g_ini
        else:
            closure_string = 'Omega_g'
            closure_variable_sym = self.sym['Omega_g']
        if closure_variable != 4:
            sub_dict[self.sym['Omega_b']] = Omega_b_ini
        else:
            closure_string = 'Omega_b'
            closure_variable_sym = self.sym['Omega_b']
        if closure_variable != 5:
            sub_dict[self.sym['Omega_c']] = Omega_c_ini
        else:
            closure_string = 'Omega_c'
            closure_variable_sym = self.sym['Omega_c']
        if closure_variable != 6:
            sub_dict[self.sym['Omega_l']] = Omega_l_ini
        else:
            closure_string = 'Omega_l'
            closure_variable_sym = self.sym['Omega_l']
        if closure_variable != 7:
            sub_dict[self.sym['Omega_n1']] = Omega_nu1_ini
        else:
            closure_string = 'Omega_n1'
            closure_variable_sym = self.sym['Omega_n1']
        if closure_variable != 8:
            sub_dict[self.sym['Omega_n2']] = Omega_nu2_ini
        else:
            closure_string = 'Omega_n2'
            closure_variable_sym = self.sym['Omega_n2']
        if closure_variable != 9:
            sub_dict[self.sym['Omega_n3']] = Omega_nu3_ini
        else:
            closure_string = 'Omega_n3'
            closure_variable_sym = self.sym['Omega_n3']
        
        for i in range(len(self.sym['K_G3_G4_syms'])):
            sub_dict[self.sym['K_G3_G4_syms'][i]] = self.params['K_G3_G4_values'][i]
        
        self.params['f_H_value'] = 1.
        sub_dict[self.sym['f_H']] = self.params['f_H_value']
        
        fried_closure = sym.simplify(self.symfunc['fried_closure'].subs(sub_dict))
        roots_raw = sym.solve(sym.Eq(fried_closure, 0), closure_variable_sym)

        if self.verbose:
            print(' -- Closure solution for %s_ini:' % closure_string, roots_raw)

        roots = []
        for root in roots_raw:
            if root.is_real:
                roots.append(root)
        
        if self.verbose:
            print(' -- Real roots for %s_ini:' % closure_string, roots)

        if bypass_closure:

            if closure_variable == 0:
                roots = [E_ini]
            elif closure_variable == 1:
                roots = [phi_ini]
            elif closure_variable == 2:
                roots = [phi_prime_ini]
            elif closure_variable == 3:
                roots = [Omega_g_ini]
            elif closure_variable == 4:
                roots = [Omega_b_ini]
            elif closure_variable == 5:
                roots = [Omega_c_ini]
            elif closure_variable == 6:
                roots = [Omega_l_ini]
            elif closure_variable == 7:
                roots = [Omega_nu1_ini]
            elif closure_variable == 8:
                roots = [Omega_nu2_ini]
            elif closure_variable == 9:
                roots = [Omega_nu3_ini]

            if self.verbose:
                print(' -- !! Bypassing closure relation solutions !!')
                print(' -- Bypassing real roots solution for %s_ini'% closure_string, roots)

        self.output['initialiser']['closure_variable'] = closure_variable
        self.output['initialiser']['closure_string'] = closure_string
        self.output['initialiser']['all_roots'] = roots_raw
        self.output['initialiser']['roots'] = roots

        if len(self.output['initialiser']['roots']) > 0:
            self.output['initialiser']['success'] = True
            self.output['success'] = True
        else:
            self.output['initialiser']['success'] = False
            self.output['success'] = False
        
        self.output['K_G3_G4_variables'] = [str(sym) for sym in self.sym['K_G3_G4_syms']],
        self.output['K_G3_G4_values'] = self.params['K_G3_G4_values']
    

    def _run_solver_ODE_NONE(self):
        """
        Returns None for solver outputs.
        """
        self.output['solver_success'] = False
        self.output['H0'] = None
        self.output['fH'] = None
        self.output['Omega_g0'] = None
        self.output['Omega_b0'] = None
        self.output['Omega_c0'] = None
        self.output['Omega_l0'] = None
        self.output['Omega_nu10'] = None
        self.output['Omega_nu20'] = None
        self.output['Omega_nu30'] = None
        self.output['Ehat'] = None
        self.output['phihat'] = None
        self.output['phihat_prime'] = None
        self.output['E'] = None
        self.output['phi'] = None
        self.output['phi_prime'] = None
        self.output['Omega_g'] = None
        self.output['Omega_b'] = None
        self.output['Omega_c'] = None
        self.output['Omega_l'] = None
        self.output['Omega_nu1'] = None
        self.output['Omega_nu2'] = None
        self.output['Omega_nu3'] = None
        self.output['w_nu1'] = None
        self.output['w_nu2'] = None
        self.output['w_nu3'] = None
    

    def _run_solver_ODE_GR(self):
        """
        Returns the GR LCDM solver outputs.
        """
        self.params['H0'] = self.params['H0_ref']
        self.params['Omega_g0'] = self.params['Omega_g0_ref']
        self.params['Omega_b0'] = self.params['Omega_b0_ref']
        self.params['Omega_c0'] = self.params['Omega_c0_ref']
        self.params['Omega_l0'] = self.params['Omega_l0_LCDM']
        self.params['Omega_nu10'] = self.params['Omega_nu10_ref']
        self.params['Omega_nu20'] = self.params['Omega_nu20_ref']
        self.params['Omega_nu30'] = self.params['Omega_nu30_ref']

        self.output['success'] = True
        self.output['solver_success'] = True
        self.output['H0'] = self.params['H0']
        self.output['fH'] = 1.
        self.output['Omega_g0'] = self.params['Omega_g0']
        self.output['Omega_b0'] = self.params['Omega_b0']
        self.output['Omega_c0'] = self.params['Omega_c0']
        self.output['Omega_l0'] = self.params['Omega_l0']
        self.output['Omega_nu10'] = self.params['Omega_nu10']
        self.output['Omega_nu20'] = self.params['Omega_nu20']
        self.output['Omega_nu30'] = self.params['Omega_nu30']

        self.output['E'] = self.compute_E_LCDM(self.output['a'], 1e-2*self.params['H0'], self.params['Omega_c0'], self.params['Omega_b0'], self.params['mnu'], w0=self.params['w0'], wa=self.params['wa'])
        self.output['phi'] = np.zeros(len(self.output['z']))
        self.output['phi_prime'] = np.zeros(len(self.output['z']))
        self.output['Omega_g'] = self.get_Omega_g(self.output['a'], 1e-2*self.params['H0'], self.output['E'])
        self.output['Omega_b'] = self.get_Omega_b(self.output['a'], self.params['Omega_b0'], self.output['E'])
        self.output['Omega_c'] = self.get_Omega_c(self.output['a'], self.params['Omega_c0'], self.output['E'])
        self.output['Omega_l'] = self.get_Omega_l(self.output['a'], self.params['Omega_l0'], self.output['E'], w0=self.params['w0'], wa=self.params['wa'])
        self.output['Omega_nu1'] = self.get_Omega_nu(self.output['a'], 1e-2*self.params['H0'], self.params['mnu'][0], self.output['E'])
        self.output['Omega_nu2'] = self.get_Omega_nu(self.output['a'], 1e-2*self.params['H0'], self.params['mnu'][1], self.output['E'])
        self.output['Omega_nu3'] = self.get_Omega_nu(self.output['a'], 1e-2*self.params['H0'], self.params['mnu'][2], self.output['E'])
        self.output['w_nu1'] = self.compute_w_nu(self.output['a'], self.params['mnu'][0])
        self.output['w_nu2'] = self.compute_w_nu(self.output['a'], self.params['mnu'][1])
        self.output['w_nu3'] = self.compute_w_nu(self.output['a'], self.params['mnu'][2])
    

    def _run_solver_ODE_HG(
            self, E_ini, phi_ini, phi_prime_ini, Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini, 
            Omega_nu1_ini, Omega_nu2_ini, Omega_nu3_ini, method, timeout):
        """
        Returns the Horndeski solver outputs.

        Parameters
        ----------
        method : str, optional
            solve_ivp method for numerical integration, use 'RK45' for general settings but switch to 'LSODA' if the solver hangs.
        timeout : float, optional
            Time in seconds to force the solver to exit and return nan, this has been added to prvent `solve_ivp` from hanging due 
            to certain variables approaching infinity. You can use different solvers, see method keyword arguement, but this will
            only work if the reason for the failure is due to the equations becoming stiff.
        """
        if self.output['success'] == False:

            self._run_solver_ODE_NONE()
        
        else:
            
            from scipy.integrate import solve_ivp

            x_arr = self.output['x']
            a_arr = self.output['a']
            z_start = self.output['z'][0]

            x_start = x_arr[0]
            x_final = x_arr[-1]

            roots = self.output['initialiser']['roots']
            closure_variable = self.output['initialiser']['closure_variable']
            
            solver_success = [False for r in roots]
            Ehat_arr = np.zeros((len(roots), len(x_arr)))
            phihat_arr = np.zeros((len(roots), len(x_arr)))
            phihat_prime_arr = np.zeros((len(roots), len(x_arr)))
            E_arr = np.zeros((len(roots), len(x_arr)))
            phi_arr = np.zeros((len(roots), len(x_arr)))
            phi_prime_arr = np.zeros((len(roots), len(x_arr)))
            Omega_g_arr = np.zeros((len(roots), len(x_arr)))
            Omega_b_arr = np.zeros((len(roots), len(x_arr)))
            Omega_c_arr = np.zeros((len(roots), len(x_arr)))
            Omega_l_arr = np.zeros((len(roots), len(x_arr)))
            Omega_nu1_arr = np.zeros((len(roots), len(x_arr)))
            Omega_nu2_arr = np.zeros((len(roots), len(x_arr)))
            Omega_nu3_arr = np.zeros((len(roots), len(x_arr)))
            w_nu1_arr = np.zeros((len(roots), len(x_arr)))
            w_nu2_arr = np.zeros((len(roots), len(x_arr)))
            w_nu3_arr = np.zeros((len(roots), len(x_arr)))
            
            H0 = np.zeros(len(roots)) 
            Omega_g0 = np.zeros(len(roots)) 
            Omega_b0 = np.zeros(len(roots)) 
            Omega_c0 = np.zeros(len(roots)) 
            Omega_l0 = np.zeros(len(roots))
            Omega_nu10 = np.zeros(len(roots)) 
            Omega_nu20 = np.zeros(len(roots)) 
            Omega_nu30 = np.zeros(len(roots)) 
            fH = np.ones(len(roots))

            for (idx, root) in enumerate(roots):
                    
                if closure_variable == 0:
                    E_ini = root
                elif closure_variable == 1:
                    phi_ini = root
                elif closure_variable == 2:
                    phi_prime_ini = root
                elif closure_variable == 3:
                    Omega_g_ini = root
                elif closure_variable == 4:
                    Omega_b_ini = root
                elif closure_variable == 5:
                    Omega_c_ini = root
                elif closure_variable == 6:
                    Omega_l_ini = root
                elif closure_variable == 7:
                    Omega_nu1_ini = root
                elif closure_variable == 8:
                    Omega_nu2_ini = root
                elif closure_variable == 9:
                    Omega_nu3_ini = root

                x_ini = x_start
                Y_ini = [E_ini, phi_ini, phi_prime_ini, Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini, Omega_nu1_ini, Omega_nu2_ini, Omega_nu3_ini]

                self._initiate_solver_status()
                self._start_timer()

                ans = solve_ivp(
                    self._compute_primes, [x_ini, x_final], Y_ini, t_eval=x_arr, method=method, 
                    args=(self.params['K_G3_G4_values'], self.params['f_H_value'], timeout), 
                    rtol = 1e-15
                )
                
                solver_success[idx] = self._solver_success

                ans = ans["y"].T
                _E_arr = ans[:,0]
                _phi_arr = ans[:,1]
                _phi_prime_arr = ans[:,2]
                _Omega_g_arr = ans[:,3]
                _Omega_b_arr = ans[:,4]
                _Omega_c_arr = ans[:,5]
                _Omega_l_arr = ans[:,6]
                _Omega_nu1_arr = ans[:,7]
                _Omega_nu2_arr = ans[:,8]
                _Omega_nu3_arr = ans[:,9]
                _w_nu1_arr = self.compute_w_nu(a_arr, self.params['mnu'][0])
                _w_nu2_arr = self.compute_w_nu(a_arr, self.params['mnu'][1])
                _w_nu3_arr = self.compute_w_nu(a_arr, self.params['mnu'][2])

                if len(_E_arr) != len(x_arr):
                    self._solver_success = False

                if self._solver_success == False:
                    split = len(_phi_prime_arr)
                    
                if self._solver_success:
                    Ehat_arr[idx] = _E_arr
                    phihat_arr[idx] = _phi_arr
                    phihat_prime_arr[idx] = _phi_prime_arr
                    Omega_g_arr[idx] = _Omega_g_arr
                    Omega_b_arr[idx] = _Omega_b_arr
                    Omega_c_arr[idx] = _Omega_c_arr
                    Omega_l_arr[idx] = _Omega_l_arr
                    Omega_nu1_arr[idx] = _Omega_nu1_arr
                    Omega_nu2_arr[idx] = _Omega_nu2_arr
                    Omega_nu3_arr[idx] = _Omega_nu3_arr
                    w_nu1_arr[idx] = _w_nu1_arr
                    w_nu2_arr[idx] = _w_nu2_arr
                    w_nu3_arr[idx] = _w_nu3_arr
                else:
                    Ehat_arr[idx][:split] = _E_arr
                    Ehat_arr[idx][split:] = np.nan
                    phihat_arr[idx][:split] = _phi_arr
                    phihat_arr[idx][split:] = np.nan
                    phihat_prime_arr[idx][:split] = _phi_prime_arr
                    phihat_prime_arr[idx][split:] = np.nan
                    Omega_g_arr[idx][:split] = _Omega_g_arr
                    Omega_g_arr[idx][split:] = np.nan
                    Omega_b_arr[idx][:split] = _Omega_b_arr
                    Omega_b_arr[idx][split:] = np.nan
                    Omega_c_arr[idx][:split] = _Omega_c_arr
                    Omega_c_arr[idx][split:] = np.nan
                    Omega_l_arr[idx][:split] = _Omega_l_arr
                    Omega_l_arr[idx][split:] = np.nan
                    Omega_nu1_arr[idx][:split] = _Omega_nu1_arr
                    Omega_nu1_arr[idx][split:] = np.nan
                    Omega_nu2_arr[idx][:split] = _Omega_nu2_arr
                    Omega_nu2_arr[idx][split:] = np.nan
                    Omega_nu3_arr[idx][:split] = _Omega_nu3_arr
                    Omega_nu3_arr[idx][split:] = np.nan
                    w_nu1_arr[idx][:split] = _w_nu1_arr
                    w_nu1_arr[idx][split:] = np.nan
                    w_nu2_arr[idx][:split] = _w_nu2_arr
                    w_nu2_arr[idx][split:] = np.nan
                    w_nu3_arr[idx][:split] = _w_nu3_arr
                    w_nu3_arr[idx][split:] = np.nan
                
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
                phi_arr[idx] = phihat_arr[idx]*self.params['f_H_value']
                phi_prime_arr[idx] = phihat_prime_arr[idx]*self.params['f_H_value']

                fH[idx] = self.params['f_H_value']
                H0[idx] = self.params['H0']

                if z_start == 0.:
                    self.params['Omega_g0'] = Omega_g_arr[idx][0]
                    self.params['Omega_b0'] = Omega_b_arr[idx][0]
                    self.params['Omega_c0'] = Omega_c_arr[idx][0]
                    self.params['Omega_l0'] = Omega_l_arr[idx][0]
                    self.params['Omega_nu10'] = Omega_nu1_arr[idx][0]
                    self.params['Omega_nu20'] = Omega_nu2_arr[idx][0]
                    self.params['Omega_nu30'] = Omega_nu3_arr[idx][0]
                else:
                    if self._solver_success == True:
                        self.params['Omega_g0'] = Omega_g_arr[idx][-1]
                        self.params['Omega_b0'] = Omega_b_arr[idx][-1]
                        self.params['Omega_c0'] = Omega_c_arr[idx][-1]
                        self.params['Omega_l0'] = Omega_l_arr[idx][-1]
                        self.params['Omega_nu10'] = Omega_nu1_arr[idx][-1]
                        self.params['Omega_nu20'] = Omega_nu2_arr[idx][-1]
                        self.params['Omega_nu30'] = Omega_nu3_arr[idx][-1]
                    else:
                        self.params['Omega_g0'] = np.nan
                        self.params['Omega_b0'] = np.nan
                        self.params['Omega_c0'] = np.nan
                        self.params['Omega_l0'] = np.nan
                        self.params['Omega_nu10'] = np.nan
                        self.params['Omega_nu20'] = np.nan
                        self.params['Omega_nu30'] = np.nan
                
                Omega_g0[idx] = self.params['Omega_g0']
                Omega_b0[idx] = self.params['Omega_b0']
                Omega_c0[idx] = self.params['Omega_c0']
                Omega_l0[idx] = self.params['Omega_l0']
                Omega_nu10[idx] = self.params['Omega_nu10']
                Omega_nu20[idx] = self.params['Omega_nu20']
                Omega_nu30[idx] = self.params['Omega_nu30']
            
            self.output['solver_success'] = solver_success
            self.output['success'] = solver_success
            self.output['H0'] = H0
            self.output['fH'] = fH
            self.output['Omega_g0'] = Omega_g0
            self.output['Omega_b0'] = Omega_b0
            self.output['Omega_c0'] = Omega_c0
            self.output['Omega_l0'] = Omega_l0
            self.output['Omega_nu10'] = Omega_nu10
            self.output['Omega_nu20'] = Omega_nu20
            self.output['Omega_nu30'] = Omega_nu30
            self.output['Ehat'] = Ehat_arr
            self.output['phihat'] = phihat_arr
            self.output['phihat_prime'] = phihat_prime_arr
            self.output['E'] = E_arr
            self.output['phi'] = phi_arr
            self.output['phi_prime'] = phi_prime_arr
            self.output['Omega_g'] = Omega_g_arr
            self.output['Omega_b'] = Omega_b_arr
            self.output['Omega_c'] = Omega_c_arr
            self.output['Omega_l'] = Omega_l_arr
            self.output['Omega_nu1'] = Omega_nu1_arr
            self.output['Omega_nu2'] = Omega_nu2_arr
            self.output['Omega_nu3'] = Omega_nu3_arr
            self.output['w_nu1'] = w_nu1_arr
            self.output['w_nu2'] = w_nu2_arr
            self.output['w_nu3'] = w_nu3_arr
    

    def _run_solver_derived_NONE(self):
        """
        Returns None for derived outputs.
        """
        self.output['H'] = None
        self.output['Dc'] = None
        self.output['G_G_4/G_N'] = None
        self.output['E_prime/E'] = None
        self.output['E_prime/E_LCDM'] = None
        self.output['E_prime'] = None
        self.output['phi_primeprime'] = None
        self.output['A'] = None
        self.output['Omega_phi'] = None
        self.output['Omega_DE'] = None
        self.output['w_DE'] = None
        self.output['Omega_phi_via_closure'] = None
        self.output['Omega_l_LCDM'] = None
        self.output['Omega_g_prime'] = None
        self.output['Omega_b_prime'] = None
        self.output['Omega_c_prime'] = None
        self.output['Omega_l_prime'] = None
        self.output['Omega_l_prime_LCDM'] = None
        self.output['Omega_nu1_prime'] = None
        self.output['Omega_nu2_prime'] = None
        self.output['Omega_nu3_prime'] = None
        self.output['w_nu1'] = None
        self.output['w_nu2'] = None
        self.output['w_nu3'] = None
        self.output['calB'] = None
        self.output['calC'] = None
        self.output['beta'] = None
        self.output['chi/delta'] = None
        self.output['M_star_sq'] = None
        self.output['alpha_M'] = None
        self.output['alpha_B'] = None
        self.output['alpha_B_prime'] = None
        self.output['alpha_K'] = None
        self.output['rho_phi'] = None
        self.output['P_phi'] = None
        self.output['w_phi'] = None
        self.output['D'] = None
        self.output['Q_s'] = None
        self.output['c_s_sq_D'] = None
        self.output['c_s_sq'] = None
        self.output['stable'] = None
    

    def _run_solver_derived_GR(self):
        """
        Returns derived outputs for LCDM GR.
        """

        x_arr = self.output['x']
        a_arr = self.output['a']
        E_arr = self.output['E']
        Omega_g_arr = self.output['Omega_g']
        Omega_b_arr = self.output['Omega_b']
        Omega_c_arr = self.output['Omega_c']
        Omega_l_arr = self.output['Omega_l']
        Omega_nu1_arr = self.output['Omega_nu1']
        Omega_nu2_arr = self.output['Omega_nu2']
        Omega_nu3_arr = self.output['Omega_nu3']

        G_G_4_G_N = np.ones(len(x_arr))

        E_prime_arr = self.compute_E_prime_LCDM(a_arr, self.params['H0']*1e-2, self.params['Omega_c0'], self.params['Omega_b0'], mnu=self.params['mnu'], w0=self.params['w0'], wa=self.params['wa'])
        
        E_prime_E_LCDM_arr = E_prime_arr/E_arr

        Omega_l_prime_arr = self.compute_Omega_l_prime(a_arr, Omega_l_arr, E_arr, E_prime_arr)

        E_prime_E_arr = E_prime_E_LCDM_arr

        phi_primeprime_arr = np.zeros(len(x_arr))

        A_arr = np.zeros(len(x_arr))

        Omega_phi_arr = np.zeros(len(x_arr))
        Omega_DE_arr = 1. - Omega_g_arr - Omega_b_arr - Omega_c_arr - Omega_nu1_arr - Omega_nu2_arr - Omega_nu3_arr

        Omega_phi_via_closure_arr = Omega_DE_arr - Omega_l_arr

        Omega_g_prime_arr = self.compute_Omega_g_prime(Omega_g_arr, E_arr, E_prime_arr)
        Omega_b_prime_arr = self.compute_Omega_b_prime(Omega_b_arr, E_arr, E_prime_arr)
        Omega_c_prime_arr = self.compute_Omega_c_prime(Omega_c_arr, E_arr, E_prime_arr)
        Omega_l_prime_arr = self.compute_Omega_l_prime(a_arr, Omega_l_arr, E_arr, E_prime_arr)
        Omega_nu1_prime_arr = self.compute_Omega_nu_prime(a_arr, Omega_nu1_arr, E_arr, E_prime_arr, self.params['H0']*1e-2, self.params['mnu'][0])
        Omega_nu2_prime_arr = self.compute_Omega_nu_prime(a_arr, Omega_nu2_arr, E_arr, E_prime_arr, self.params['H0']*1e-2, self.params['mnu'][1])
        Omega_nu3_prime_arr = self.compute_Omega_nu_prime(a_arr, Omega_nu3_arr, E_arr, E_prime_arr, self.params['H0']*1e-2, self.params['mnu'][2])

        calB_arr = np.zeros(len(x_arr))
        calC_arr = np.zeros(len(x_arr))
        beta_arr = np.zeros(len(x_arr))

        chioverdelta_arr = np.zeros(len(x_arr))

        M_star_sq_arr = np.ones(len(x_arr))*(self.params['mass_ratios']['M_p']**2.)
        alpha_M_arr = np.zeros(len(x_arr))
        alpha_B_arr = np.zeros(len(x_arr))
        alpha_B_prime_arr = np.zeros(len(x_arr))
        alpha_K_arr = np.zeros(len(x_arr))
        
        rho_phi_arr = np.zeros(len(x_arr))
        P_phi_arr = np.zeros(len(x_arr))
        w_phi_arr = np.zeros(len(x_arr))

        w_DE_arr = self.compute_w_l(a_arr, w0=self.params['w0'], wa=self.params['wa'])*Omega_l_arr + w_phi_arr*Omega_phi_arr
        w_DE_arr /= Omega_l_arr + Omega_phi_arr
        
        D_arr = np.zeros(len(x_arr))
        Q_s_arr = np.zeros(len(x_arr))
        c_s_sq_D_arr = np.nan*np.ones(len(x_arr))
        c_s_sq_arr = np.nan*np.ones(len(x_arr))
        stable = True

        H_arr = self.output['H0']*E_arr

        from scipy.integrate import cumulative_trapezoid
        
        f = 1./E_arr
        if self.output['initialiser']['forwards']:
            Dc_arr = self.const['Dh']*cumulative_trapezoid(f[::-1], x=self.output['z'][::-1], initial=0.)[::-1]
        else:
            Dc_arr = self.const['Dh']*cumulative_trapezoid(f, x=self.output['z'], initial=0.)

        self.output['H'] = H_arr
        self.output['Dc'] = Dc_arr
        self.output['G_G_4/G_N'] = G_G_4_G_N
        self.output['E_prime/E'] = E_prime_E_arr
        self.output['E_prime/E_LCDM'] = E_prime_E_LCDM_arr
        self.output['E_prime'] = E_prime_arr
        self.output['phi_primeprime'] = phi_primeprime_arr
        self.output['A'] = A_arr
        self.output['Omega_phi'] = Omega_phi_arr
        self.output['Omega_DE'] = Omega_DE_arr
        self.output['w_DE'] = w_DE_arr
        self.output['Omega_phi_via_closure'] = Omega_phi_via_closure_arr 
        self.output['Omega_l_LCDM'] = Omega_l_arr
        self.output['Omega_g_prime'] = Omega_g_prime_arr
        self.output['Omega_b_prime'] = Omega_b_prime_arr
        self.output['Omega_c_prime'] = Omega_c_prime_arr
        self.output['Omega_l_prime'] = Omega_l_prime_arr
        self.output['Omega_nu1_prime'] = Omega_nu1_prime_arr
        self.output['Omega_nu2_prime'] = Omega_nu2_prime_arr
        self.output['Omega_nu3_prime'] = Omega_nu3_prime_arr
        self.output['calB'] = calB_arr
        self.output['calC'] = calC_arr
        self.output['beta'] = beta_arr
        self.output['chi/delta'] = chioverdelta_arr
        self.output['M_star_sq'] = M_star_sq_arr
        self.output['alpha_M'] = alpha_M_arr
        self.output['alpha_B'] = alpha_B_arr
        self.output['alpha_B_prime'] = alpha_B_prime_arr
        self.output['alpha_K'] = alpha_K_arr
        self.output['rho_phi'] = rho_phi_arr
        self.output['P_phi'] = P_phi_arr
        self.output['w_phi'] = w_phi_arr
        self.output['D'] = D_arr
        self.output['Q_s'] = Q_s_arr
        self.output['c_s_sq_D'] = c_s_sq_D_arr
        self.output['c_s_sq'] = c_s_sq_arr
        self.output['stable'] = stable


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
            phi_arr = self.output['phi']
            phi_prime_arr = self.output['phi_prime']
            Omega_g_arr = self.output['Omega_g']
            Omega_b_arr = self.output['Omega_b']
            Omega_c_arr = self.output['Omega_c']
            Omega_l_arr = self.output['Omega_l']
            Omega_nu1_arr = self.output['Omega_nu1']
            Omega_nu2_arr = self.output['Omega_nu2']
            Omega_nu3_arr = self.output['Omega_nu3']
            w_nu1_arr = self.output['w_nu1']
            w_nu2_arr = self.output['w_nu2']
            w_nu3_arr = self.output['w_nu3']

            roots = self.output['initialiser']['roots']
            
            G_G_4_G_N = np.zeros((len(roots), len(x_arr)))
            E_prime_E_arr = np.zeros((len(roots), len(x_arr)))
            E_prime_arr = np.zeros((len(roots), len(x_arr)))
            phi_primeprime_arr = np.zeros((len(roots), len(x_arr)))
            A_arr = np.zeros((len(roots), len(x_arr)))
            Omega_phi_arr = np.zeros((len(roots), len(x_arr)))
            Omega_DE_arr = np.zeros((len(roots), len(x_arr)))
            w_DE_arr = np.zeros((len(roots), len(x_arr)))
            Omega_phi_via_closure_arr = np.zeros((len(roots), len(x_arr)))
            Omega_g_prime_arr = np.zeros((len(roots), len(x_arr)))
            Omega_b_prime_arr = np.zeros((len(roots), len(x_arr)))
            Omega_c_prime_arr = np.zeros((len(roots), len(x_arr)))
            Omega_l_prime_arr = np.zeros((len(roots), len(x_arr)))
            Omega_nu1_prime_arr = np.zeros((len(roots), len(x_arr)))
            Omega_nu2_prime_arr = np.zeros((len(roots), len(x_arr)))
            Omega_nu3_prime_arr = np.zeros((len(roots), len(x_arr)))
            calB_arr = np.zeros((len(roots), len(x_arr)))
            calC_arr = np.zeros((len(roots), len(x_arr)))
            beta_arr = np.zeros((len(roots), len(x_arr)))
            chioverdelta_arr = np.zeros((len(roots), len(x_arr)))
            M_star_sq_arr = np.zeros((len(roots), len(x_arr)))
            alpha_M_arr = np.zeros((len(roots), len(x_arr)))
            alpha_B_arr = np.zeros((len(roots), len(x_arr)))
            alpha_B_prime_arr = np.zeros((len(roots), len(x_arr)))
            alpha_K_arr = np.zeros((len(roots), len(x_arr)))
            rho_phi_arr = np.zeros((len(roots), len(x_arr)))
            P_phi_arr = np.zeros((len(roots), len(x_arr)))
            w_phi_arr = np.zeros((len(roots), len(x_arr)))
            D_arr = np.zeros((len(roots), len(x_arr)))
            Q_s_arr = np.zeros((len(roots), len(x_arr)))
            c_s_sq_D_arr = np.zeros((len(roots), len(x_arr)))
            c_s_sq_arr = np.zeros((len(roots), len(x_arr)))

            stable = [True for r in roots]

            for (idx, _) in enumerate(roots):

                self.params['H0'] = self.output['H0'][idx]
                self.params['f_H_value'] = self.output['fH'][idx]

                variables = [
                    E_arr[idx], phi_arr[idx], phi_prime_arr[idx], Omega_g_arr[idx], Omega_b_arr[idx], Omega_c_arr[idx], Omega_l_arr[idx], 
                    Omega_nu1_arr[idx], w_nu1_arr[idx], Omega_nu2_arr[idx], w_nu2_arr[idx], Omega_nu3_arr[idx], w_nu3_arr[idx], 
                    self.compute_w_l(a_arr, w0=self.params['w0'], wa=self.params['wa']), *self.params['K_G3_G4_values'], self.params['f_H_value']]

                G_G_4_G_N[idx] = self.lambda_funcs['G_G_4/G_N'](*variables)
                E_prime_E_arr[idx] = self.lambda_funcs['E_prime/E_lambda'](*variables)
                E_prime_arr[idx] = E_prime_E_arr[idx] * E_arr[idx]

                variables = [E_prime_arr[idx], *variables]

                phi_primeprime_arr[idx] = self.lambda_funcs['phi_primeprime_lambda'](*variables)

                variables = [phi_primeprime_arr[idx], *variables]

                A_arr[idx] = self.lambda_funcs['A_lambda'](*variables)
                Omega_phi_arr[idx] = self.lambda_funcs['Omega_phi_lambda'](*variables)
                Omega_DE_arr[idx] = 1. - Omega_g_arr[idx] - Omega_b_arr[idx] - Omega_c_arr[idx] - Omega_nu1_arr[idx] - Omega_nu2_arr[idx] - Omega_nu3_arr[idx]
                Omega_phi_via_closure_arr[idx] = 1 - Omega_g_arr[idx] - Omega_b_arr[idx] - Omega_c_arr[idx] - Omega_l_arr[idx] - Omega_nu1_arr[idx] - Omega_nu2_arr[idx] - Omega_nu3_arr[idx]

                Omega_g_prime_arr[idx] = self.compute_Omega_g_prime(Omega_g_arr[idx], E_arr[idx], E_prime_arr[idx])
                Omega_b_prime_arr[idx] = self.compute_Omega_b_prime(Omega_b_arr[idx], E_arr[idx], E_prime_arr[idx])
                Omega_c_prime_arr[idx] = self.compute_Omega_c_prime(Omega_c_arr[idx], E_arr[idx], E_prime_arr[idx])
                Omega_l_prime_arr[idx] = self.compute_Omega_l_prime(a_arr, Omega_l_arr[idx], E_arr[idx], E_prime_arr[idx])
                Omega_nu1_prime_arr[idx] = self.compute_Omega_nu_prime(a_arr, Omega_nu1_arr[idx], E_arr[idx], E_prime_arr[idx], self.params['H0']*1e-2, self.params['mnu'][0])
                Omega_nu2_prime_arr[idx] = self.compute_Omega_nu_prime(a_arr, Omega_nu2_arr[idx], E_arr[idx], E_prime_arr[idx], self.params['H0']*1e-2, self.params['mnu'][1])
                Omega_nu3_prime_arr[idx] = self.compute_Omega_nu_prime(a_arr, Omega_nu3_arr[idx], E_arr[idx], E_prime_arr[idx], self.params['H0']*1e-2, self.params['mnu'][2])

                calB_arr[idx] = self.lambda_funcs['calB_lambda'](*variables)
                calC_arr[idx] = self.lambda_funcs['calC_lambda'](*variables)
                beta_arr[idx] = self.lambda_funcs['beta_lambda'](*variables)
                chioverdelta_arr[idx] = self.compute_chi_over_delta(a_arr, E_arr[idx], calB_arr[idx], calC_arr[idx])

                M_star_sq_arr[idx] = self.lambda_funcs['M_star_sq'](*variables)
                alpha_M_arr[idx] = self.lambda_funcs['alpha_M'](*variables)
                alpha_B_arr[idx] = self.lambda_funcs['alpha_B'](*variables)
                alpha_K_arr[idx] = self.lambda_funcs['alpha_K'](*variables)
                
                if self.output['initialiser']['forwards']:
                    alpha_B_prime_arr[idx] = np.gradient(alpha_B_arr[idx], x_arr)
                else:
                    alpha_B_prime_arr[idx] = np.gradient(alpha_B_arr[idx][::-1], x_arr[::-1])[::-1]
                
                rho_phi_arr[idx] = self.lambda_funcs['rho_phi'](self.params['H0'], *variables)
                P_phi_arr[idx] = self.lambda_funcs['P_phi'](self.params['H0'], *variables)
                w_phi_arr[idx] = P_phi_arr[idx]/rho_phi_arr[idx]
                
                w_DE_arr[idx] = self.compute_w_l(a_arr, w0=self.params['w0'], wa=self.params['wa'])*Omega_l_arr[idx] + w_phi_arr[idx]*Omega_phi_arr[idx]
                w_DE_arr[idx] /= Omega_l_arr[idx] + Omega_phi_arr[idx]

                D_arr[idx] = self.lambda_funcs['D'](*variables)
                Q_s_arr[idx] = self.lambda_funcs['Q_s'](*variables)
                c_s_sq_D_arr[idx] = self.lambda_funcs['c_s_sq_D'](alpha_B_prime_arr[idx], *variables)
                c_s_sq_arr[idx] = self.lambda_funcs['c_s_sq'](alpha_B_prime_arr[idx], *variables)

                if Q_s_arr.all() > 0 and c_s_sq_arr.all() > 0:
                    stable[idx] = True
                else:
                    stable[idx] = False

            H_arr = np.zeros_like(E_arr)
            Dc_arr = np.zeros_like(E_arr)
            
            for (idx, _) in enumerate(roots):
                
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
            self.output['E_prime/E'] = E_prime_E_arr
            self.output['E'] = E_arr
            self.output['E_prime'] = E_prime_arr
            self.output['phi_primeprime'] = phi_primeprime_arr
            self.output['A'] = A_arr
            self.output['Omega_phi'] = Omega_phi_arr
            self.output['Omega_DE'] = Omega_DE_arr
            self.output['w_DE'] = w_DE_arr
            self.output['Omega_phi_via_closure'] = Omega_phi_via_closure_arr 
            self.output['Omega_g_prime'] = Omega_g_prime_arr
            self.output['Omega_b_prime'] = Omega_b_prime_arr
            self.output['Omega_c_prime'] = Omega_c_prime_arr
            self.output['Omega_l_prime'] = Omega_l_prime_arr
            self.output['Omega_nu1_prime'] = Omega_nu1_prime_arr
            self.output['Omega_nu2_prime'] = Omega_nu2_prime_arr
            self.output['Omega_nu3_prime'] = Omega_nu3_prime_arr
            self.output['calB'] = calB_arr
            self.output['calC'] = calC_arr
            self.output['beta'] = beta_arr
            self.output['chi/delta'] = chioverdelta_arr
            self.output['M_star_sq'] = M_star_sq_arr
            self.output['alpha_M'] = alpha_M_arr
            self.output['alpha_B'] = alpha_B_arr
            self.output['alpha_B_prime'] = alpha_B_prime_arr
            self.output['alpha_K'] = alpha_K_arr
            self.output['rho_phi'] = rho_phi_arr
            self.output['P_phi'] = P_phi_arr
            self.output['w_phi'] = w_phi_arr
            self.output['D'] = D_arr
            self.output['Q_s'] = Q_s_arr
            self.output['c_s_sq_D'] = c_s_sq_D_arr
            self.output['c_s_sq'] = c_s_sq_arr
            self.output['stable'] = stable
    

    def _reverse_outputs(self, derived):
        """
        Reverse order of outputted quantities, so that functions start from the early universe to late.
        """
        
        if derived:
            keys = [
                'a', 'z', 'x', 'Ehat', 'phihat', 'phihat_prime', 'E', 'phi', 'phi_prime', 
                'Omega_g', 'Omega_b', 'Omega_c', 'Omega_l', 'Omega_nu1', 'Omega_nu2', 'Omega_nu3', 
                'w_nu1', 'w_nu2', 'w_nu3', 
                # derived quantities
                'H', 'Dc', 'G_G_4/G_N', 'E_prime/E',
                'E_prime', 'phi_primeprime', 'A', 'Omega_phi', 'Omega_DE', 'w_DE', 'Omega_phi_via_closure', 
                'Omega_g_prime', 'Omega_b_prime', 'Omega_c_prime', 'Omega_l_prime', 
                'Omega_nu1_prime', 'Omega_nu2_prime', 'Omega_nu3_prime', 
                'calB', 'calC', 'beta', 'chi/delta', 'M_star_sq',
                'alpha_M', 'alpha_B', 'alpha_B_prime', 'alpha_K' ,
                'rho_phi', 'P_phi', 'w_phi', 'D', 'Q_s', 'c_s_sq_D', 'c_s_sq'
            ]
        else:
            keys = [
                'a', 'z', 'x', 'Ehat', 'phihat', 'phihat_prime', 'E', 'phi', 'phi_prime', 
                'Omega_g', 'Omega_b', 'Omega_c', 'Omega_l', 'Omega_nu1', 'Omega_nu2', 'Omega_nu3', 
                'w_nu1', 'w_nu2', 'w_nu3'
            ]

        for key in keys:
            if self.output[key] is not None:
                if self.output[key].ndim == 1:
                    self.output[key] = self.output[key][::-1]
                elif self.output[key].ndim == 2:
                    self.output[key] = self.output[key][:,::-1]


    def run_solver(
            self, z_max=1100., Npoints=1000, forwards=True, GR=False, closure_variable=2, 
            phi_ini=1e-5, phi_prime_ini=0.9, method='RK45', timeout=5, derived=True, LCDM_ini=True, 
            values_ini=None, bypass_closure=False
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
        closure_variable : str, optional
            Variable used to set the initial conditions.
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
            Directly supply initial values for initial conditions for [E, Omega_g, Omega_b, Omega_c, Omega_l, Omega_nu1, Omega_nu2, Omega_nu3]
        bypass_closure : bool, optional
            Bypass closure solver to set quantities directly, this means the closure relation will not be satisfied, use with care!
        
        Returns
        -------
        outputs : dict
            Dictionary containing numerical solver solutions.
        """

        if self.verbose:
            print('Hi-COLA: Running numerical ODE solver')
            print(' - Initialising solver...')

        self.output = {}
        E_ini, Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini, Omega_nu1_ini, Omega_nu2_ini, Omega_nu3_ini = self._run_solver_start(z_max, Npoints, forwards)

        if LCDM_ini == False:
            assert values_ini is not None, "values_ini must be a list with the initial values for [E_ini, Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini, Omega_nu1, Omega_nu2, Omega_nu3]"
            E_ini, Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini, Omega_nu1_ini, Omega_nu2_ini, Omega_nu3_ini = values_ini[0], values_ini[1], values_ini[2], values_ini[3], values_ini[4], values_ini[5], values_ini[6], values_ini[7]

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

            self._run_solver_HG_root_finder(
                closure_variable, E_ini, phi_ini, phi_prime_ini, Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini, 
                Omega_nu1_ini, Omega_nu2_ini, Omega_nu3_ini, bypass_closure)

            if self.output['success']:
                
                if self.verbose:
                    print(' - Numerically solving ODEs for expansion history and scalar field evolution.')

                self._run_solver_ODE_HG(
                    E_ini, phi_ini, phi_prime_ini, Omega_g_ini, Omega_b_ini, Omega_c_ini, Omega_l_ini, 
                    Omega_nu1_ini, Omega_nu2_ini, Omega_nu3_ini, method, timeout)

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

            # Compute equation of state via alternative equation # TODO: remove? -- this for the moment only provides a sanity check
            self.compute_w_phi()

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
        
        if np.isscalar(self.output['Omega_c0']):

            fname_expansion = fname_prefix + 'expansion.txt'
            data = np.column_stack([self.output['a'], self.output['E'], self.output['E_prime/E']])
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
                data = np.column_stack([self.output['a'], self.output['E'][idx], self.output['E_prime/E'][idx]])
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
