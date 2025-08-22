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

    def __init__(self):
        """
        Initialises the Horndeski model class.
        """
        # Enable print statements
        self.verbose = True
        # Construct dictionary containing sympy symbolic variables.
        self.sym = {}
        self.sym['a'] = sym.symbols('a')
        self.sym['E'] = sym.symbols('E')
        self.sym['Eprime'] = sym.symbols('Eprime')
        self.sym['phi'] = sym.symbols('phi')
        self.sym['phiprime'] = sym.symbols('phiprime')
        self.sym['phiprimeprime'] = sym.symbols('phiprimeprime')
        self.sym['X'] = sym.symbols('X')
        self.sym['M_pG4'] = sym.symbols('M_{pG4}')
        self.sym['M_KG4'] = sym.symbols('M_{KG4}')
        self.sym['M_G3s'] = sym.symbols('M_{G3s}')
        self.sym['M_sG4'] = sym.symbols('M_{sG4}')
        self.sym['M_G3G4'] = sym.symbols('M_{G3G4}')
        self.sym['M_Ks'] = sym.symbols('M_{Ks}')
        self.sym['M_gp'] = sym.symbols('M_{gp}')
        self.sym['Omega_r'] = sym.symbols('Omega_r')
        self.sym['Omega_m'] = sym.symbols('Omega_m')
        self.sym['Omega_l'] = sym.symbols('Omega_l')
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
        self.params['mass_ratios'] = {
            'M_pG4': 1, 
            'M_KG4': 1, 
            'M_G3s': 1, 
            'M_sG4': 1, 
            'M_G3G4': 1, 
            'M_Ks': 1, 
            'M_gp': 1,
            'M_sp': 1,
        }
        self.outputs = {}
        self.timer = {'t0': None}
    

    def _start_timer(self):
        self.timer['t0'] = time.time()
    

    def _check_timer(self):
        timenow = time.time()
        return timenow - self.timer['t0']
    

    def clean_params(self):
        """
        Reinitialised the class parameters.
        """
        self.params = {}
        self.outputs = {}
    

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


    def get_Omega_phi(self):
        """
        Returns the scalar field fractional density function, following equation 2.4 in https://arxiv.org/pdf/2209.01666.pdf.
        """
        # Note: replaces the omega_phi function
        term1 = (self.sym['Omega_m'] + self.sym['Omega_r'] + self.sym['Omega_l'])*((self.sym['M_pG4']**2.)/(2.*self.symfunc['G4']) - 1.)
        term2 = (self.sym['M_KG4']**2.)*self.sym['X']*self.symfunc['Kx']/(self.sym['E']**2.) - (self.sym['M_KG4']**2.)*self.symfunc['K']/(2.*(self.sym['E']**2.))
        term2 += 3*self.sym['M_G3G4']*self.sym['M_sG4']*self.sym['X']*self.sym['phiprime']*self.symfunc['G3x']
        term2 += -self.sym['M_G3G4']*self.sym['M_sG4']*self.sym['X']*self.symfunc['G3phi']/(self.sym['E']**2.) - 3*self.sym['phiprime']*self.symfunc['G4phi']
        self.symfunc['Omega_phi'] = term1 + (1/(3.*self.symfunc['G4']))*term2


    def get_EprimeE(self):
        """
        Returns the evolution of the dimensionless Hubble function E=H/H0, as a function of log(a), follows equation 2.5 
        https://arxiv.org/pdf/2209.01666.pdf.
        """
        # Note: replaces EprimeEODERHS function.
        M_G4s = 1./self.sym['M_sG4']
        A = (self.sym['M_Ks']**2.)*self.symfunc['Kx'] - self.sym['M_G3s']*self.symfunc['G3phi'] + 2.*self.sym['X']*self.sym['M_G3s']*self.symfunc['G3phix']
        A += 6*(self.sym['E']**2)*self.sym['phiprime']*(self.sym['M_G3s']*self.symfunc['G3x'] + self.sym['X']*self.sym['M_G3s']*self.symfunc['G3xx']) 
        A += (self.sym['E']**2.)*(self.sym['phiprime']**2.)*((self.sym['M_Ks']**2.)*self.symfunc['Kxx'] - 2.*self.sym['M_G3s']*self.symfunc['G3phix'])

        B1 = 6.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3x'] - 6.*(M_G4s**2)*self.symfunc['G4phi']

        B2 = 3.*self.sym['phiprime']*((self.sym['M_Ks']**2.)*self.symfunc['Kx'] - 2.*self.sym['M_G3s']*self.symfunc['G3phi'] + 2.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3phix'])
        B2 += (self.sym['phiprime']**2.)*((self.sym['M_Ks']**2.)*self.symfunc['Kxphi'] - 2.*self.sym['M_G3s']*self.symfunc['G3phiphi']) - ((self.sym['M_Ks']**2.)/(self.sym['E']**2.))*self.symfunc['Kphi']
        B2 += -12.*(M_G4s**2.)*self.symfunc['G4phi'] + 18.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3x'] + 2.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3phiphi']*(1./self.sym['E']**2.)

        term1 = 1.+ (1./2.)*self.sym['M_G3G4']*self.sym['M_sG4']*self.sym['X']*self.symfunc['G3x']*(1./self.symfunc['G4'])*(B1/A) - (self.symfunc['G4phi']*B1)/(2*self.symfunc['G4']*A)

        term21 = (self.sym['M_KG4']**2.)*self.symfunc['K']*(1./(self.sym['E']**2.)) - 2.*self.sym['M_sG4']*self.sym['M_G3G4']*(1./(self.sym['E']**2.))*self.sym['X']*self.symfunc['G3phi']
        term22 = 4.*self.symfunc['G4phi']*self.sym['phiprime'] + 4.*self.sym['X']*self.symfunc['G4phiphi']*(1./(self.sym['E']**2.))
        term2 = (-1./(4.*self.symfunc['G4']))*(term21 + term22)

        term3 = (-1./2.)*(((self.sym['Omega_r'] - 3.*self.sym['Omega_l'])/(2.*self.symfunc['G4']))*(self.sym['M_pG4']**2.) + 3.) 
        term4 = (-1./2.)*self.sym['M_G3G4']*self.sym['M_sG4']*self.sym['X']*self.symfunc['G3x']*(B2/(self.symfunc['G4']*A)) + (self.symfunc['G4phi']*B2)/(2*self.symfunc['G4']*A)
        RHS = term2 + term3 + term4

        self.symfunc['EprimeE'] =  sym.simplify((RHS)/sym.simplify(term1))


    def get_EprimeE_safe(self):
        """
        This is identical to get_EprimeE, except the variable "A" is replaced by a 
        constant, "threshold". This function is used if a given model has a tendency to send
        A close to 0, and computational errors take it over 0. By replacing A with "threshold"
        this behaviour is (artificially) avoided.
        """
        M_G4s = 1./self.sym['M_sG4']
        A = self.sym['threshold']*self.sym['threshold_sign']

        B1 = 6.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3x'] - 6.*(M_G4s**2)*self.symfunc['G4phi']

        B2 = 3.*self.sym['phiprime']*((self.sym['M_Ks']**2.)*self.symfunc['Kx'] - 2.*self.sym['M_G3s']*self.symfunc['G3phi'] + 2.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3phix'])
        B2 += (self.sym['phiprime']**2.)*((self.sym['M_Ks']**2.)*self.symfunc['Kxphi'] - 2.*self.sym['M_G3s']*self.symfunc['G3phiphi']) - ((self.sym['M_Ks']**2.)/(self.sym['E']**2.))*self.symfunc['Kphi']
        B2 += -12.*(M_G4s**2.)*self.symfunc['G4phi'] + 18.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3x'] + 2.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3phiphi']*(1./self.sym['E']**2.)

        term1 = 1.+ (1./2.)*self.sym['M_G3G4']*self.sym['M_sG4']*self.sym['X']*self.symfunc['G3x']*(1./self.symfunc['G4'])*(B1/A) - (self.symfunc['G4phi']*B1)/(2*self.symfunc['G4']*A)

        term21 = (self.sym['M_KG4']**2.)*self.symfunc['K']*(1./(self.sym['E']**2.)) - 2.*self.sym['M_sG4']*self.sym['M_G3G4']*(1./(self.sym['E']**2.))*self.sym['X']*self.symfunc['G3phi']
        term22 = 4.*self.symfunc['G4phi']*self.sym['phiprime'] + 4.*self.sym['X']*self.symfunc['G4phiphi']*(1./(self.sym['E']**2.))
        term2 = (-1./(4.*self.symfunc['G4']))*(term21 + term22)

        term3 = (-1./2.)*((self.sym['Omega_r']/(2.*self.symfunc['G4']))*(self.sym['M_pG4']**2.) - 3.*self.sym['Omega_l'] + 3. )
        term4 = (-1./2.)*self.sym['M_G3G4']*self.sym['M_sG4']*self.sym['X']*self.symfunc['G3x']*(B2/(self.symfunc['G4']*A)) + (self.symfunc['G4phi']*B2)/(2*self.symfunc['G4']*A)
        RHS = term2 + term3 + term4

        self.symfunc['EprimeE_safe'] =  sym.simplify((RHS)/sym.simplify(term1))

    
    def get_phiprimeprime(self):
        """
        Return the phiprimeprime function, equation 2.6 in https://arxiv.org/abs/2209.01666
        """
        # Note: replaces the phiprimeprimeODERHS function
        M_G4s = 1./self.sym['M_sG4']
        A = (self.sym['M_Ks']**2.)*self.symfunc['Kx'] - self.sym['M_G3s']*self.symfunc['G3phi'] + 2.*self.sym['X']*self.sym['M_G3s']*self.symfunc['G3phix']
        A += 6*(self.sym['E']**2)*self.sym['phiprime']*(self.sym['M_G3s']*self.symfunc['G3x'] + self.sym['X']*self.sym['M_G3s']*self.symfunc['G3xx']) + (self.sym['E']**2.)*(self.sym['phiprime']**2.)*((self.sym['M_Ks']**2.)*self.symfunc['Kxx'] - 2.*self.sym['M_G3s']*self.symfunc['G3phix'])
        B1 = 6.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3x'] - 6.*(M_G4s**2)*self.symfunc['G4phi']
        B2 = 3.*self.sym['phiprime']*((self.sym['M_Ks']**2.)*self.symfunc['Kx'] - 2.*self.sym['M_G3s']*self.symfunc['G3phi'] + 2.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3phix'])
        B2 += (self.sym['phiprime']**2.)*((self.sym['M_Ks']**2.)*self.symfunc['Kxphi'] - 2.*self.sym['M_G3s']*self.symfunc['G3phiphi']) - ((self.sym['M_Ks']**2.)/(self.sym['E']**2.))*self.symfunc['Kphi']
        B2 += -12.*(M_G4s**2.)*self.symfunc['G4phi'] + 18.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3x'] + 2.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3phiphi']*(1./(self.sym['E']**2.))
        B =  B1*(self.sym['Eprime']/self.sym['E'])+B2
        self.symfunc['phiprimeprime'] = -1.*((B/A) + (self.sym['Eprime']/self.sym['E'])*self.sym['phiprime'])
        

    def get_phiprimeprime_safe(self):
        """
        This is identical to get_phiprimeprime function and analogous to the get_EprimeE_safe function.
        """
        # Note: replaces the phiprimeprimeODERHS function
        M_G4s = 1./self.sym['M_sG4']
        A = self.sym['threshold']*self.sym['threshold_sign']
        B1 = 6.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3x'] - 6.*(M_G4s**2)*self.symfunc['G4phi']
        B2 = 3.*self.sym['phiprime']*((self.sym['M_Ks']**2.)*self.symfunc['Kx'] - 2.*self.sym['M_G3s']*self.symfunc['G3phi'] + 2.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3phix'])
        B2 += (self.sym['phiprime']**2.)*((self.sym['M_Ks']**2.)*self.symfunc['Kxphi'] - 2.*self.sym['M_G3s']*self.symfunc['G3phiphi']) - ((self.sym['M_Ks']**2.)/(self.sym['E']**2.))*self.symfunc['Kphi']
        B2 += -12.*(M_G4s**2.)*self.symfunc['G4phi'] + 18.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3x'] + 2.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3phiphi']*(1./self.sym['E']**2.)
        B =  B1*(self.sym['Eprime']/self.sym['E'])+B2
        self.symfunc['phiprimeprime_safe'] = -1.*( (B/A)  + (self.sym['Eprime']/self.sym['E'])*self.sym['phiprime'])
        

    def get_fried_closure(self):
        """
        Returns the RHS of the Friedmann closure relation -> equation 2.3 in https://arxiv.org/abs/2209.01666, 
        should be equal to 1. 
        """
        # Note: replaces the fried_closure function
        self.get_Omega_phi()
        self.symfunc['fried_closure'] = self.symfunc['Omega_phi'] + self.sym['Omega_m'] + self.sym['Omega_r'] + self.sym['Omega_l'] - 1.
        # equation A.3 in https://arxiv.org/abs/2209.01666, 
        Xreal = (1./2.)*(self.sym['E']**2.)*(self.sym['phiprime']**2.)
        self.symfunc['fried_closure'] = self.symfunc['fried_closure'].subs(self.sym['X'], Xreal)

    
    def get_A(self):
        """
        Code equivalent of equation 2.7 in https://arxiv.org/abs/2209.01666.
        """
        # Note: replaces the A_func
        self.symfunc['A'] = (self.sym['M_Ks']**2.)*self.symfunc['Kx'] - self.sym['M_G3s']*self.symfunc['G3phi'] + 2.*self.sym['X']*self.sym['M_G3s']*self.symfunc['G3phix']
        self.symfunc['A'] += 6*(self.sym['E']**2)*self.sym['phiprime']*(self.sym['M_G3s']*self.symfunc['G3x'] + self.sym['X']*self.sym['M_G3s']*self.symfunc['G3xx']) 
        self.symfunc['A'] += (self.sym['E']**2.)*(self.sym['phiprime']**2.)*((self.sym['M_Ks']**2.)*self.symfunc['Kxx'] - 2.*self.sym['M_G3s']*self.symfunc['G3phix'])
    

    def get_B2(self):
        """
        Code equivalent of equation 2.10 in https://arxiv.org/abs/2209.01666.
        """
        # Note: replaces the B2_func
        M_G4s = 1./self.sym['M_sG4']
        self.symfunc['B2'] = 3.*self.sym['phiprime']*((self.sym['M_Ks']**2.)*self.symfunc['Kx'] - 2.*self.sym['M_G3s']*self.symfunc['G3phi'] + 2.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3phix'])
        self.symfunc['B2'] += (self.sym['phiprime']**2.)*((self.sym['M_Ks']**2.)*self.symfunc['Kxphi'] - 2.*self.sym['M_G3s']*self.symfunc['G3phiphi']) - ((self.sym['M_Ks']**2.)/(self.sym['E']**2.))*self.symfunc['Kphi']
        self.symfunc['B2'] += -12.*(M_G4s**2.)*self.symfunc['G4phi'] + 18.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3x'] + 2.*self.sym['M_G3s']*self.sym['X']*self.symfunc['G3phiphi']*(1./self.sym['E']**2.)


    def get_theta(self):
        """
        Code equivalent of equation 2.15 in https://arxiv.org/abs/2209.01666
        """
        # Note: replaces theta function
        # All terms seemed to be multiplied by E, not so in eq. 2.15, why?
        term1 = self.sym['M_sG4']*self.sym['M_G3G4']*self.sym['E']*self.sym['phiprime']*self.sym['X']*self.symfunc['G3x']/self.sym['M_pG4']/self.sym['M_pG4']
        term2 = 2.*self.sym['E']*self.symfunc['G4']/self.sym['M_pG4']/self.sym['M_pG4']
        term3 = self.sym['E']*self.sym['phiprime']*self.symfunc['G4phi']/self.sym['M_pG4']/self.sym['M_pG4']
        self.symfunc['theta'] = -1.*term1 + term2 + term3
        Xreal = 0.5*(self.sym['E']**2.)*self.sym['phiprime']**2.
        self.symfunc['theta'] = self.symfunc['theta'].subs(self.sym['X'], Xreal)
        

    def get_calE(self):
        """
        Code equivalent of equation 5 in https://arxiv.org/abs/1111.6749.
        """
        # Note: replaces calE function
        term1 = 2.*self.sym['M_KG4']*self.sym['M_KG4']*self.sym['X']*self.symfunc['Kx']/self.sym['M_pG4']/self.sym['M_pG4']
        term2 = self.sym['M_KG4']*self.sym['M_KG4']*self.symfunc['K']/self.sym['M_pG4']/self.sym['M_pG4']
        term3 = 6.*self.sym['M_sG4']*self.sym['M_G3G4']*self.sym['E']*self.sym['E']*self.sym['X']*self.symfunc['G3x']*self.sym['phiprime']/self.sym['M_pG4']/self.sym['M_pG4']
        term4 = 2.*self.sym['M_sG4']*self.sym['M_G3G4']*self.sym['X']*self.symfunc['G3phi']/self.sym['M_pG4']/self.sym['M_pG4']
        term5 = 6.*self.sym['E']*self.sym['E']*self.symfunc['G4']/self.sym['M_pG4']/self.sym['M_pG4']
        term6 = 6.*self.sym['E']*self.sym['E']*self.symfunc['G4phi']*self.sym['phiprime']/self.sym['M_pG4']/self.sym['M_pG4']
        self.symfunc['calE'] = term1 - term2 + term3 - term4 - term5 - term6


    def get_calP(self):
        """
        Code equivalent of equation 6 in https://arxiv.org/abs/1111.6749
        """
        # Note: replaces calP function
        term1 = self.sym['M_KG4']*self.sym['M_KG4']*self.symfunc['K']/self.sym['M_pG4']/self.sym['M_pG4']
        term2coeff = 2.*self.sym['M_sG4']*self.sym['X']/self.sym['M_pG4']
        term2bracket = self.sym['M_G3s']*self.symfunc['G3phi'] + self.sym['E']*self.sym['M_G3s']*self.symfunc['G3x']*(self.sym['Eprime']*self.sym['phiprime'] + self.sym['E']*self.sym['phiprimeprime'])
        term2 = term2coeff*term2bracket
        term3 = 2.*self.symfunc['G4']*(3.*self.sym['E']*self.sym['E'] + 2*self.sym['E']*self.sym['Eprime'])/self.sym['M_pG4']/self.sym['M_pG4']
        term4coeff = 2*self.symfunc['G4phi']/self.sym['M_pG4']/self.sym['M_pG4']
        term4bracket = self.sym['E']*(self.sym['Eprime']*self.sym['phiprime'] + self.sym['E']*self.sym['phiprimeprime']) + 2*self.sym['E']*self.sym['E']*self.sym['phiprime']
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
        self.sym['dtheta_dphiprime'] = sym.diff(self.symfunc['theta'], self.sym['phiprime'])
        self.sym['thetaprime'] = self.sym['dtheta_dE']*self.sym['Eprime'] + self.sym['dtheta_dphiprime']*self.sym['phiprimeprime']

        A0 = self.sym['thetaprime']/self.sym['E'] + self.symfunc['theta']/self.sym['E'] 
        A0 += -2.*self.symfunc['G4']/self.sym['M_pG4']/self.sym['M_pG4'] - 4.*self.symfunc['G4phi']*self.sym['phiprime']/self.sym['M_pG4']/self.sym['M_pG4'] 
        A0 += -(self.symfunc['calE'] + self.symfunc['calP'])/(2.*self.sym['E']*self.sym['E'])
        self.symfunc['alpha0'] = self.sym['M_pG4']*self.sym['M_pG4']*A0/2./self.symfunc['G4']
        

    def get_alpha1(self):
        """
        Code equivalent of equation 3.5 (see 2.16-2.19) in https://arxiv.org/abs/2209.01666.
        Given by alpha1 = A1/2G4.
        """
        # Note: replaces the alpha1 function
        A1 = 2.*self.symfunc['G4phi']*self.sym['phiprime']/self.sym['M_pG4']/self.sym['M_pG4']
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
        B0 = self.sym['M_sG4']*self.sym['M_G3G4']*self.sym['X']*self.symfunc['G3x']*self.sym['phiprime']/self.sym['M_pG4']/self.sym['M_pG4']
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

    
    def set_mass_ratios(self, M_pG4=1., M_KG4=1., M_G3s=1., M_sG4=1., M_G3G4=1., M_Ks=1., M_gp=1.):
        """
        Allows the users to assign specific values to mass ratios, those unchanged will be set to 1.

        Parameters
        ----------
        M_pG4 : float, optional
            TODO
        M_KG4 : float, optional
            TODO
        M_G3s : float, optional
            TODO
        M_sG4 : float, optional
            TODO
        M_G3G4 : float, optional
            TODO
        M_Ks : float, optional
            TODO
        M_gp : float, optional
            TODO
        """
        self.params['mass_ratios'] = {
            'M_pG4': M_pG4, 
            'M_KG4': M_KG4, 
            'M_G3s': M_G3s, 
            'M_sG4': M_sG4, 
            'M_G3G4': M_G3G4, 
            'M_Ks': M_Ks, 
            'M_gp': M_gp ,
            'M_sp': M_sG4/M_pG4,
        }


    def construct_model(self):
        """
        Constructs Horndeski model with user defined functions.
        """

        if self._check_symfunc_keys(['K', 'G3', 'G4']) == False:
            assert False, 'Functions for K, G3 and G4 remain undefined.'
        else:
            self._get_K_G3_G4_syms()
            self.get_K_derivatives()
            self.get_G3_derivatives()
            self.get_G4_derivatives()

            Xreal = 0.5*(self.sym['E']**2.)*self.sym['phiprime']**2.

            sub_dict = {
                self.sym['X']: Xreal,
                self.sym['M_pG4']: self.params['mass_ratios']['M_pG4'],
                self.sym['M_KG4']: self.params['mass_ratios']['M_KG4'],
                self.sym['M_G3s']: self.params['mass_ratios']['M_G3s'],
                self.sym['M_sG4']: self.params['mass_ratios']['M_sG4'],
                self.sym['M_G3G4']: self.params['mass_ratios']['M_G3G4'],
                self.sym['M_Ks']: self.params['mass_ratios']['M_Ks'],
                self.sym['M_gp']: self.params['mass_ratios']['M_gp']
            }

            self.get_EprimeE()
            EprimeE = self.symfunc['EprimeE'].subs(sub_dict)

            self.get_EprimeE_safe()
            EprimeE_safe = self.symfunc['EprimeE_safe'].subs(sub_dict)

            self.get_phiprimeprime()
            phiprimeprime = self.symfunc['phiprimeprime'].subs(sub_dict)

            self.get_phiprimeprime_safe()
            phiprimeprime_safe = self.symfunc['phiprimeprime_safe'].subs(sub_dict)

            self.get_A()
            A = self.symfunc['A'].subs(sub_dict)

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

            # we will copy these substituted and simplified functions to the class symfunc dictionary, we avoided 
            # doing this before as some of these are re-called and redefined in the 'get' functions.
            self.symfunc['EprimeE'] = EprimeE
            self.symfunc['EprimeE_safe'] = EprimeE_safe
            self.symfunc['phiprimeprime'] = phiprimeprime
            self.symfunc['phiprimeprime_safe'] = phiprimeprime_safe
            self.symfunc['A'] = A
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

            # Lambdify functions
            self.lambda_funcs['B2_lambda'] = sym.lambdify([self.sym['E'], self.sym['phiprime'], *self.sym['K_G3_G4_syms']], self.symfunc['B2'], "scipy")
            self.lambda_funcs['fried_closure_lambda'] = sym.lambdify(
                [self.sym['E'], self.sym['phiprime'], self.sym['Omega_r'], self.sym['Omega_m'], self.sym['Omega_l'], *self.sym['K_G3_G4_syms']],
                self.symfunc['fried_closure']
            )
            self.lambda_funcs['EprimeE_lambda'] = sym.lambdify([self.sym['E'], self.sym['phiprime'], self.sym['Omega_r'], self.sym['Omega_l'], *self.sym['K_G3_G4_syms']], self.symfunc['EprimeE'], "scipy")
            self.lambda_funcs['EprimeE_safe_lambda'] = sym.lambdify(
                [self.sym['E'], self.sym['phiprime'], self.sym['Omega_r'], self.sym['Omega_l'], self.sym['threshold'], self.sym['threshold_sign'], *self.sym['K_G3_G4_syms']],
                self.symfunc['EprimeE_safe'], "scipy"
            )
            self.lambda_funcs['phiprimeprime_lambda'] = sym.lambdify([self.sym['E'], self.sym['Eprime'], self.sym['phiprime'], *self.sym['K_G3_G4_syms']], self.symfunc['phiprimeprime'], "scipy")
            self.lambda_funcs['phiprimeprime_safe_lambda'] = sym.lambdify(
                [self.sym['E'], self.sym['Eprime'], self.sym['phiprime'], self.sym['threshold'], self.sym['threshold_sign'], *self.sym['K_G3_G4_syms']], 
                self.symfunc['phiprimeprime_safe'], "scipy"
            )

            self.lambda_funcs['Omega_phi_lambda'] = sym.lambdify([self.sym['E'], self.sym['phiprime'], *self.sym['K_G3_G4_syms']], self.symfunc['Omega_phi'])
            self.lambda_funcs['A_lambda'] = sym.lambdify([self.sym['E'], self.sym['phiprime'], *self.sym['K_G3_G4_syms']], self.symfunc['A'], "scipy")
            self.lambda_funcs['alpha0_lambda'] = sym.lambdify([self.sym['E'], self.sym['Eprime'], self.sym['phiprime'], self.sym['phiprimeprime'], *self.sym['K_G3_G4_syms']], self.symfunc['alpha0'])
            self.lambda_funcs['alpha1_lambda'] = sym.lambdify([self.sym['E'], self.sym['Eprime'], self.sym['phiprime'], self.sym['phiprimeprime'], *self.sym['K_G3_G4_syms']], self.symfunc['alpha1'])
            self.lambda_funcs['alpha2_lambda'] = sym.lambdify([self.sym['E'], self.sym['Eprime'], self.sym['phiprime'], self.sym['phiprimeprime'], *self.sym['K_G3_G4_syms']], self.symfunc['alpha2'])
            self.lambda_funcs['beta0_lambda'] = sym.lambdify([self.sym['E'], self.sym['Eprime'], self.sym['phiprime'], self.sym['phiprimeprime'], *self.sym['K_G3_G4_syms']], self.symfunc['beta0'])
            self.lambda_funcs['calB_lambda'] = sym.lambdify([self.sym['E'], self.sym['Eprime'], self.sym['phiprime'], self.sym['phiprimeprime'], *self.sym['K_G3_G4_syms']], self.symfunc['calB'])
            self.lambda_funcs['calC_lambda'] = sym.lambdify([self.sym['E'], self.sym['Eprime'], self.sym['phiprime'], self.sym['phiprimeprime'], *self.sym['K_G3_G4_syms']], self.symfunc['calC'])
            self.lambda_funcs['beta_lambda'] = sym.lambdify([self.sym['E'], self.sym['Eprime'], self.sym['phiprime'], self.sym['phiprimeprime'], *self.sym['K_G3_G4_syms']], self.symfunc['beta'])
    

    def set_cosmo_params(self, H0_ref, Omega_m0_ref, Omega_r0_ref, fphi, K_G3_G4_values):
        """
        Set cosmological and Horndeski parameters.

        Parameters
        ----------
        H0_ref : float
            Reference LCDM Hubble constant.
        Omega_m0_ref : float
            Reference LCDM matter density at redshift zero.
        Omega_r0_ref : float
            Reference LCDM radiation density at redshift zero.
        fphi : float
            The fraction of the scalar field density as a fraction of the full dark energy density (including a cosmological constant).
        K_G3_G4_values : list
            A list of values for the Horndeski specific variables. This must match the length of the user defined variable. 
            Check self.sym['K_G3_G4_syms'] to see what variables are expected.
        """
        self.params['H0_ref'] = H0_ref
        self.params['Omega_m0_ref'] = Omega_m0_ref
        self.params['Omega_r0_ref'] = Omega_r0_ref
        self.params['fphi'] = fphi
        self.params['Omega_phi0_ref'] = fphi*(1. - self.params['Omega_r0_ref'] - self.params['Omega_m0_ref'])
        self.params['Omega_l0_ref'] = 1. - self.params['Omega_r0_ref'] - self.params['Omega_m0_ref'] - self.params['Omega_phi0_ref']
        self.params['Omega_l0_LCDM'] = 1. - self.params['Omega_r0_ref'] - self.params['Omega_m0_ref']
        assert len(K_G3_G4_values) == len(self.sym['K_G3_G4_syms']), "Length of Horndeski K_G3_G4_values must match number of defined K, G3, G4 variables."
        self.params['K_G3_G4_values'] = K_G3_G4_values
    
    
    def _fried_closure_wrapper(self, cl_val, cl_var, fried_closure_lambda, E_ini, phi_prime_ini, Omega_r_ini, Omega_m_ini, Omega_l_ini, K_G3_G4_values):
        """
        Wrapper function for the Friedmann closure relation.

        Parameters
        ----------
        cl_val : float
            Closure variable value.
        cl_var : int
            Closure variable index, 0 = E, 1 = phi_prime, 2 = Omega_r, 3 = Omega_m, 4 = Omega_l
        fried_closure_lamba : func
            Function for the closure relation, this should just `self.lambda_funcs['fried_closure_lambda']`.
        E_ini : float
            E initial value, if cl_var = 0 this value is computed through the closure relation.
        phi_prime_ini : float
            phi_prime initial value, if cl_var = 1 this value is computed through the closure relation.
        Omega_r_ini : float
            Omega radiation initial value, if cl_var = 2 this value is computed through the closure relation.
        Omega_m_prime_ini : float
            Omega matter initial value, if cl_var = 3 this value is computed through the closure relation.
        Omega_l_ini : float
            Omega lambda initial value, if cl_var = 4 this value is computed through the closure relation.
        K_G3_G4_values : float
            A list of values for the Horndeski specific variables. This must match the length of the user defined variable. 
            Check self.sym['K_G3_G4_syms'] to see what variables are expected.
        """
        assert cl_var >= 0 and cl_var <= 4, "Closure variable unsupported, must be between 0 and 4 inclusive."
        if cl_var == 0:
            #Closure used to compute E0
            return fried_closure_lambda(cl_val, phi_prime_ini, Omega_r_ini, Omega_m_ini, Omega_l_ini, *K_G3_G4_values)
        if cl_var == 1:
            #Closure used to compute phi0
            return fried_closure_lambda(E_ini, cl_val, Omega_r_ini, Omega_m_ini, Omega_l_ini, *K_G3_G4_values)
        if cl_var == 2:
            #Closure used to compute Omega_r0
            return fried_closure_lambda(E_ini, phi_prime_ini, cl_val, Omega_m_ini, Omega_l_ini, *K_G3_G4_values)
        if cl_var == 3:
            #Closure used to compute Omega_m0
            return fried_closure_lambda(E_ini, phi_prime_ini, Omega_r_ini, cl_val, Omega_l_ini, *K_G3_G4_values)
        if cl_var == 4:
            #Closure used to compute Omega_l0
            return fried_closure_lambda(E_ini, phi_prime_ini, Omega_r_ini, Omega_m_ini, cl_val, *K_G3_G4_values)
    

    def compute_Omega_r_prime(self, Omega_r, E, E_prime):
        """
        Computes Omega radiation prime.

        Parameters
        ----------
        Omega_r : float or array
            Radiation density.
        E : float or array
            Normalised Hubble expansion.
        E_prime : float or array
            Derivative of the normalised Hubble expansion.
        
        Returns
        -------
        Omega_r_prime : float or array
            Omega radiation prime.
        """
        E_prime_E = E_prime/E
        Omega_r_prime = -Omega_r*(4. + 2.*E_prime_E)
        return Omega_r_prime


    def compute_Omega_m_prime(self, Omega_m, E, E_prime):
        """
        Computes Omega radiation prime.

        Parameters
        ----------
        Omega_m : float or array
            Matter density.
        E : float or array
            Normalised Hubble expansion.
        E_prime : float or array
            Derivative of the normalised Hubble expansion.
        
        Returns
        -------
        Omega_r_prime : float or array
            Omega matter prime.
        """
        E_prime_E = E_prime/E
        Omega_m_prime = -Omega_m*(3. + 2.*E_prime_E)
        return Omega_m_prime


    def compute_Omega_l_prime(self, Omega_l, E, E_prime):
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
        """
        E_prime_E = E_prime/E
        Omega_l_prime = -2.*Omega_l*E_prime_E
        return Omega_l_prime


    def compute_chi_over_delta(self, a, E, calB, calC):
        """
        Computes the chi/delta function.

        Parameters
        ----------
        a : float or array
            Scale factor.
        E : float or array
            Normalised Hubble expansion.
        calB : float or array
            TODO
        calC : float or array
            TODO

        Return
        ------
        chioverdelta : float or array   
            TODO
        """
        chioverdelta = calB * calC * self.params['Omega_m0']/((E**2)*(a**3)) # TODO: there's a G_G4/G_N in 3.14 in https://arxiv.org/pdf/2209.01666 which is not included here...
        return chioverdelta
    
    
    def _initiate_solver_status(self):
        """
        Initiate solver status monitor.
        """
        self._solver_success = None


    def _compute_primes(self, x, Y, K_G3_G4_values, timeout=5):
        """
        Compute prime functions for numerical solver.

        Parameters
        ----------
        x : float
            Current value of log(a).
        Y : list
            List containing current [phi_prime, E, Omega_r, Omega_m, Omega_l] values.
        K_G3_G4_values : list
            A list of values for the Horndeski specific variables. This must match the length of the user defined variable. 
            Check self.sym['K_G3_G4_syms'] to see what variables are expected.
        threshold : float, optional
            Numerical solver threshold to switch to 'safe' functions.
        timeout : float, optional
            Time in seconds to force the solver to fail.
        """
        # Note2self: x is unused -- solve_ivp probably requires x so cannot remove this

        # `_` used to denote current value.
        _phi_prime, _E, _Omega_r, _Omega_m, _Omega_l = Y

        E_prime_E = self.lambda_funcs['EprimeE_lambda'](_E, _phi_prime, _Omega_r, _Omega_l, *K_G3_G4_values)
        E_prime = E_prime_E*_E
        phi_primeprime = self.lambda_funcs['phiprimeprime_lambda'](_E, E_prime, _phi_prime, *K_G3_G4_values)
       
        Omega_r_prime = self.compute_Omega_r_prime(_Omega_r, _E, E_prime)
        Omega_m_prime = self.compute_Omega_m_prime(_Omega_m, _E, E_prime)
        Omega_l_prime = self.compute_Omega_l_prime(_Omega_l, _E, E_prime)

        timenow = self._check_timer()

        if timenow >= timeout:

            phi_primeprime = np.nan
            E_prime = np.nan
            Omega_r_prime = np.nan
            Omega_m_prime = np.nan
            Omega_l_prime = np.nan

        if np.isfinite([phi_primeprime, E_prime, Omega_r_prime, Omega_m_prime, Omega_l_prime]).all() == False:
            self._solver_success = False

        return [phi_primeprime, E_prime, Omega_r_prime, Omega_m_prime, Omega_l_prime]


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
        Bx = 2 + self._linear_growth_int['interp_Eprime_vs_a'](a)/self._linear_growth_int['interp_E_vs_a'](a)
        Cx = 3.*Omega_m0*mu/(2.*(self._linear_growth_int['interp_E_vs_a'](a)**2)*(a**3))
        d2D = Cx*D - Bx*D
        return [dD, d2D]
    

    def get_linear_growth(self):
        """
        Compute the first order linear growth function.
        """
        
        from scipy.interpolate import interp1d
        from scipy.integrate import solve_ivp

        # Run sanity checks to test whether linear growth functions can be computed

        if np.isscalar(self.output['Omega_m0']):

            check = True
            if self.output['E'] is None or np.isfinite(self.output['E']).all() == False:
                check = False
            if self.output['E_prime'] is None or np.isfinite(self.output['E_prime']).all() == False:
                check = False
            if self.output['beta'] is None or np.isfinite(self.output['beta']).all() == False:
                check = False

            if check:

                self._linear_growth_int = {}
                self._linear_growth_int['interp_E_vs_a'] = interp1d(self.output['a'], self.output['E'], kind='cubic', fill_value='extrapolate')
                self._linear_growth_int['interp_Eprime_vs_a'] = interp1d(self.output['a'], self.output['E_prime'], kind='cubic', fill_value='extrapolate')
                self._linear_growth_int['interp_beta_vs_a'] = interp1d(self.output['a'], self.output['beta'], kind='cubic', fill_value='extrapolate')

                # set up initial conditons, assuming matter domination.
                
                x_ini = self.output['x'][0]
                D_ini = self.output['a'][0]
                dD_ini = self.output['a'][0]
                y_ini = [D_ini, dD_ini]

                # End position
                x_final = self.output['x'][-1]
                
                # Solve forward
                ans = solve_ivp(self._linear_growth, (x_ini, x_final), y_ini, t_eval=self.output['x'], args=(self.output['Omega_m0'],))
                    
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

        else:
            
            D1 = np.zeros(np.shape(self.output['E']))
            f1 = np.zeros(np.shape(self.output['E']))

            for idx in range(0, len(self.output['Omega_m0'])):

                check = True
                if self.output['E'] is None or np.isfinite(self.output['E'][idx]).all() == False:
                    check = False
                if self.output['E_prime'] is None or np.isfinite(self.output['E_prime'][idx]).all() == False:
                    check = False
                if self.output['beta'] is None or np.isfinite(self.output['beta'][idx]).all() == False:
                    check = False

                if check:

                    self._linear_growth_int = {}
                    self._linear_growth_int['interp_E_vs_a'] = interp1d(self.output['a'], self.output['E'][idx], kind='cubic', fill_value='extrapolate')
                    self._linear_growth_int['interp_Eprime_vs_a'] = interp1d(self.output['a'], self.output['E_prime'][idx], kind='cubic', fill_value='extrapolate')
                    self._linear_growth_int['interp_beta_vs_a'] = interp1d(self.output['a'], self.output['beta'][idx], kind='cubic', fill_value='extrapolate')

                    # set up initial conditons, assuming matter domination.
                    
                    x_ini = self.output['x'][0]
                    D_ini = self.output['a'][0]
                    dD_ini = self.output['a'][0]
                    y_ini = [D_ini, dD_ini]

                    # End position
                    x_final = self.output['x'][-1]
                    
                    # Solve forward
                    ans = solve_ivp(self._linear_growth, (x_ini, x_final), y_ini, t_eval=self.output['x'], args=(self.output['Omega_m0'][idx],))
                        
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
        Bx = 2 + self._linear_growth_int['interp_Eprime_vs_a'](a)/self._linear_growth_int['interp_E_vs_a'](a)
        Cx = 3.*Omega_m0*mu/(2.*(self._linear_growth_int['interp_E_vs_a'](a)**2)*(a**3))
        d2D2 = Cx*(D2 - D1**2) - Bx*dD2
        return [dD2, d2D2]


    def get_linear_growth_2(self):
        """
        Compute the second order linear growth function assuming mu1 = mu2.
        """
        
        from scipy.interpolate import interp1d
        from scipy.integrate import solve_ivp

        # Run sanity checks to test whether linear growth functions can be computed

        if np.isscalar(self.output['Omega_m0']):

            check = True
            if self.output['E'] is None or np.isfinite(self.output['E']).all() == False:
                check = False
            if self.output['E_prime'] is None or np.isfinite(self.output['E_prime']).all() == False:
                check = False
            if self.output['beta'] is None or np.isfinite(self.output['beta']).all() == False:
                check = False
            if self.output['D1'] is None or np.isfinite(self.output['D1']).all() == False:
                check = False

            if check:
                    
                self._linear_growth_int = {}
                self._linear_growth_int['interp_E_vs_a'] = interp1d(self.output['a'], self.output['E'], kind='cubic', fill_value='extrapolate')
                self._linear_growth_int['interp_Eprime_vs_a'] = interp1d(self.output['a'], self.output['E_prime'], kind='cubic', fill_value='extrapolate')
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
                ans = solve_ivp(self._linear_growth_2, (x_ini, x_final), y_ini, t_eval=self.output['x'], args=(self.output['Omega_m0'],))
                    
                # Combine solutions
                _D2 = ans.y[0]
                _dD2 = ans.y[1]

                _f2 = _dD2/_D2

                D2, f2 = _D2, _f2
            
            else:

                D2, f2 = None, None

        else:

            D2 = np.zeros(np.shape(self.output['E']))
            f2 = np.zeros(np.shape(self.output['E']))

            for idx in range(0, len(self.output['Omega_m0'])):

                check = True
                if self.output['E'] is None or np.isfinite(self.output['E'][idx]).all() == False:
                    check = False
                if self.output['E_prime'] is None or np.isfinite(self.output['E_prime'][idx]).all() == False:
                    check = False
                if self.output['beta'] is None or np.isfinite(self.output['beta'][idx]).all() == False:
                    check = False
                if self.output['D1'] is None or np.isfinite(self.output['D1'][idx]).all() == False:
                    check = False

                if check:

                    self._linear_growth_int = {}
                    self._linear_growth_int['interp_E_vs_a'] = interp1d(self.output['a'], self.output['E'][idx], kind='cubic', fill_value='extrapolate')
                    self._linear_growth_int['interp_Eprime_vs_a'] = interp1d(self.output['a'], self.output['E_prime'][idx], kind='cubic', fill_value='extrapolate')
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
                    ans = solve_ivp(self._linear_growth_2, (x_ini, x_final), y_ini, t_eval=self.output['x'], args=(self.output['Omega_m0'][idx],))
                        
                    # Combine solutions
                    _D2 = ans.y[0]
                    _dD2 = ans.y[1]

                    _f2 = _dD2/_D2

                    D2[idx], f2[idx] = _D2, _f2
                
                else:
                    
                    D2[idx] = np.nan * np.ones(len(self.output['x']))
                    f2[idx] = np.nan * np.ones(len(self.output['x']))

        self.output['D2'] = D2
        self.output['f2'] = f2


    def run_solver(self, z_max=1000., Npoints=1000, forwards=True, GR=False, closure_variable=1, phi_prime_ini=0.9, method='RK45', 
        timeout=5, compute_growth=True):
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
        phi_prime_ini : float, optional
            phi_prime initial value, if the closure_variable=1 then this is used as an initial guess.
        threshold : float, optional
            Numerical solver threshold to switch to 'safe' functions.
        method : str, optional
            solve_ivp method for numerical integration, use 'RK45' for general settings but switch to 'LSODA' if the solver hangs.
        timeout : float, optional
            Time in seconds to force the solver to exit and return nan, this has been added to prvent `solve_ivp` from hanging due 
            to certain variables approaching infinity. You can use different solvers, see method keyword arguement, but this will
            only work if the reason for the failure is due to the equations becoming stiff.
        compute_growth : bool, optional
            Flag to instruct the solver to compute growth functions.
        
        Returns
        -------
        outputs : dict
            Dictionary containing numerical solver solutions.
        """

        # from scipy.optimize import fsolve # TODO remove
        from scipy.integrate import solve_ivp

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

        x_start = x_arr[0]
        x_final = x_arr[-1]

        z_start = z_arr[0]

        # Let's guess the values of the variables by assuming the solution lies close to the reference LCDM values.

        E_ini = lcdm.compute_Ez_LCDM(z_start, self.params['Omega_r0_ref'], self.params['Omega_m0_ref'])
        Omega_r_ini = lcdm.compute_Omega_r_z_LCDM(z_start, self.params['Omega_r0_ref'], self.params['Omega_m0_ref'])
        Omega_m_ini = lcdm.compute_Omega_m_z_LCDM(z_start, self.params['Omega_r0_ref'], self.params['Omega_m0_ref'])
        Omega_l_ini = (1.-self.params['fphi'])*lcdm.compute_Omega_l_z_LCDM(z_start, self.params['Omega_r0_ref'], self.params['Omega_m0_ref'])
        
        if GR == False:
            
            sub_dict = {}
            
            assert closure_variable >= 0 and closure_variable <= 4, "Closure variable unsupported, must be between 0 and 4 inclusive."

            if closure_variable != 0:
                sub_dict[self.sym['E']] = E_ini
            else:
                closure_variable_sym = self.sym['E']
            if closure_variable != 1:
                sub_dict[self.sym['phiprime']] = phi_prime_ini
            else:
                closure_variable_sym = self.sym['phiprime']
            if closure_variable != 2:
                sub_dict[self.sym['Omega_r']] = Omega_r_ini
            else:
                closure_variable_sym = self.sym['Omega_r']
            if closure_variable != 3:
                sub_dict[self.sym['Omega_m']] = Omega_m_ini
            else:
                closure_variable_sym = self.sym['Omega_m']
            if closure_variable != 4:
                sub_dict[self.sym['Omega_l']] = Omega_l_ini
            else:
                closure_variable_sym = self.sym['Omega_l']
            
            for i in range(len(self.sym['K_G3_G4_syms'])):
                sub_dict[self.sym['K_G3_G4_syms'][i]] = self.params['K_G3_G4_values'][i]
            
            fried_closure = sym.simplify(self.symfunc['fried_closure'].subs(sub_dict))
            
            roots_raw = sym.solve(sym.Eq(fried_closure, 0), closure_variable_sym)

            roots = []
            for root in roots_raw:
                if root.is_real:
                    roots.append(root)

            if len(roots) > 0:
                
                solver_success = [False for r in roots]
                phi_prime_arr = np.zeros((len(roots), len(x_arr)))
                E_arr = np.zeros((len(roots), len(x_arr)))
                Omega_r_arr = np.zeros((len(roots), len(x_arr)))
                Omega_m_arr = np.zeros((len(roots), len(x_arr)))
                Omega_l_arr = np.zeros((len(roots), len(x_arr)))

                H0 = [False for r in roots] 
                Omega_r0 = [False for r in roots] 
                Omega_m0 = [False for r in roots]
                Omega_l0 = [False for r in roots]

                E_prime_E_arr = np.zeros((len(roots), len(x_arr)))
                E_prime_arr = np.zeros((len(roots), len(x_arr)))
                phi_primeprime_arr = np.zeros((len(roots), len(x_arr)))
                A_arr = np.zeros((len(roots), len(x_arr)))
                Omega_phi_arr = np.zeros((len(roots), len(x_arr)))
                Omega_DE_arr = np.zeros((len(roots), len(x_arr)))
                Omega_phi_via_closure_arr = np.zeros((len(roots), len(x_arr)))
                Omega_r_prime_arr = np.zeros((len(roots), len(x_arr)))
                Omega_m_prime_arr = np.zeros((len(roots), len(x_arr)))
                Omega_l_prime_arr = np.zeros((len(roots), len(x_arr)))
                calB_arr = np.zeros((len(roots), len(x_arr)))
                calC_arr = np.zeros((len(roots), len(x_arr)))
                beta_arr = np.zeros((len(roots), len(x_arr)))
                chioverdelta_arr = np.zeros((len(roots), len(x_arr)))

                # Not sure the LCDM arrays are necessary...
                E_prime_E_LCDM_arr = lcdm.compute_EprimeE_x_LCDM(x_arr, self.params['Omega_r0_ref'], self.params['Omega_m0_ref'])
                Omega_l_LCDM_arr = lcdm.compute_Omega_l_x_LCDM(x_arr, self.params['Omega_r0_ref'], self.params['Omega_m0_ref'])
                Omega_l_prime_LCDM_arr = lcdm.compute_Omega_l_prime_LCDM(E_prime_E_LCDM_arr, 1-self.params['Omega_r0_ref']-self.params['Omega_m0_ref'])
            
                for (idx, root) in enumerate(roots):
                    
                    if closure_variable == 0:
                        E_ini = root
                    elif closure_variable == 1:
                        phi_prime_ini = root
                    elif closure_variable == 2:
                        Omega_r_ini = root
                    elif closure_variable == 3:
                        Omega_m_ini = root
                    elif closure_variable == 4:
                        Omega_l_ini = root

                    x_ini = x_start
                    Y_ini = [phi_prime_ini, E_ini, Omega_r_ini, Omega_m_ini, Omega_l_ini]

                    self._initiate_solver_status()
                    self._start_timer()

                    ans = solve_ivp(
                        self._compute_primes, [x_ini, x_final], Y_ini, t_eval=x_arr, method=method, 
                        args=(self.params['K_G3_G4_values'], timeout), 
                        rtol = 1e-15
                    )

                    if self._solver_success is None:
                        self._solver_success = True
                    
                    solver_success[idx] = self._solver_success

                    ans = ans["y"].T
                    _phi_prime_arr = ans[:,0]
                    _E_arr = ans[:,1]
                    _Omega_r_arr = ans[:,2]
                    _Omega_m_arr = ans[:,3]
                    _Omega_l_arr = ans[:,4]

                    if len(_phi_prime_arr) != len(x_arr):
                        self._solver_success = False
                    
                    if self._solver_success == False:
                        split = len(_phi_prime_arr)
                    
                    if self._solver_success:
                        phi_prime_arr[idx] = _phi_prime_arr
                        E_arr[idx] = _E_arr
                        Omega_r_arr[idx] = _Omega_r_arr
                        Omega_m_arr[idx] = _Omega_m_arr
                        Omega_l_arr[idx] = _Omega_l_arr
                    else:
                        phi_prime_arr[idx][:split] = _phi_prime_arr
                        phi_prime_arr[idx][split:] = np.nan
                        E_arr[idx][:split] = _E_arr
                        E_arr[idx][split:] = np.nan
                        Omega_r_arr[idx][:split] = _Omega_r_arr
                        Omega_r_arr[idx][split:] = np.nan
                        Omega_m_arr[idx][:split] = _Omega_m_arr
                        Omega_m_arr[idx][split:] = np.nan
                        Omega_l_arr[idx][:split] = _Omega_l_arr
                        Omega_l_arr[idx][split:] = np.nan

                    if z_start == 0.:
                        self.params['H0'] = self.params['H0_ref']*E_arr[idx][0]
                        self.params['Omega_r0'] = Omega_r_arr[idx][0]
                        self.params['Omega_m0'] = Omega_m_arr[idx][0]
                        self.params['Omega_l0'] = Omega_l_arr[idx][0]
                    else:
                        if self._solver_success == True:
                            self.params['H0'] = self.params['H0_ref']*E_arr[idx][-1]
                            self.params['Omega_r0'] = Omega_r_arr[idx][-1]
                            self.params['Omega_m0'] = Omega_m_arr[idx][-1]
                            self.params['Omega_l0'] = Omega_l_arr[idx][-1]
                        else:
                            self.params['H0'] = np.nan
                            self.params['Omega_r0'] = np.nan
                            self.params['Omega_m0'] = np.nan
                            self.params['Omega_l0'] = np.nan
                    
                    H0[idx] = self.params['H0']
                    Omega_r0[idx] = self.params['Omega_r0']
                    Omega_m0[idx] = self.params['Omega_m0']
                    Omega_l0[idx] = self.params['Omega_l0']

                    E_prime_E_arr[idx] = self.lambda_funcs['EprimeE_lambda'](E_arr[idx], phi_prime_arr[idx], Omega_r_arr[idx], Omega_l_arr[idx], *self.params['K_G3_G4_values'])
                    E_prime_arr[idx] = E_prime_E_arr[idx] * E_arr[idx]
                
                    phi_primeprime_arr[idx] = self.lambda_funcs['phiprimeprime_lambda'](E_arr[idx], E_prime_arr[idx], phi_prime_arr[idx], *self.params['K_G3_G4_values'])

                    A_arr[idx] = self.lambda_funcs['A_lambda'](E_arr[idx], phi_prime_arr[idx], *self.params['K_G3_G4_values'])

                    Omega_phi_arr[idx] = self.lambda_funcs['Omega_phi_lambda'](E_arr[idx], phi_prime_arr[idx], *self.params['K_G3_G4_values'])

                    Omega_DE_arr[idx] = 1. - Omega_r_arr[idx] - Omega_m_arr[idx]

                    Omega_phi_via_closure_arr[idx] = 1 - Omega_r_arr[idx] - Omega_m_arr[idx] - Omega_l_arr[idx]

                    Omega_r_prime_arr[idx] = self.compute_Omega_r_prime(Omega_r_arr[idx], E_arr[idx], E_prime_arr[idx])
                    Omega_m_prime_arr[idx] = self.compute_Omega_m_prime(Omega_m_arr[idx], E_arr[idx], E_prime_arr[idx])
                    Omega_l_prime_arr[idx] = self.compute_Omega_l_prime(Omega_l_arr[idx], E_arr[idx], E_prime_arr[idx])

                    calB_arr[idx] = self.lambda_funcs['calB_lambda'](E_arr[idx], E_prime_arr[idx], phi_prime_arr[idx], phi_primeprime_arr[idx], *self.params['K_G3_G4_values'])
                    calC_arr[idx] = self.lambda_funcs['calC_lambda'](E_arr[idx], E_prime_arr[idx], phi_prime_arr[idx], phi_primeprime_arr[idx], *self.params['K_G3_G4_values'])
                    beta_arr[idx] = self.lambda_funcs['beta_lambda'](E_arr[idx], E_prime_arr[idx], phi_prime_arr[idx], phi_primeprime_arr[idx], *self.params['K_G3_G4_values'])
                    chioverdelta_arr[idx] = self.compute_chi_over_delta(a_arr, E_arr[idx], calB_arr[idx], calC_arr[idx])

            else:
                roots = None
                solver_success = None
                phi_prime_arr = None
                E_arr = None
                Omega_r_arr = None
                Omega_m_arr = None
                Omega_l_arr = None
                x_arr = None
                a_arr = None
                z_arr = None
                H0 = None
                Omega_r0 = None
                Omega_m0 = None
                Omega_l0 = None
                E_prime_E_LCDM_arr = None
                Omega_l_LCDM_arr = None
                Omega_l_prime_LCDM_arr = None
                E_prime_E_arr = None
                E_prime_arr = None
                phi_primeprime_arr = None
                A_arr = None
                Omega_phi_arr = None
                Omega_DE_arr = None
                Omega_phi_via_closure_arr = None
                Omega_r_prime_arr = None
                Omega_m_prime_arr = None
                Omega_l_prime_arr = None
                calB_arr = None
                calC_arr = None
                beta_arr = None
                chioverdelta_arr = None
        else:
            
            self.params['H0'] = self.params['H0_ref']
            H0 = self.params['H0']
            self.params['Omega_r0'] = self.params['Omega_r0_ref']
            Omega_r0 = self.params['Omega_r0']
            self.params['Omega_m0'] = self.params['Omega_m0_ref']
            Omega_m0 = self.params['Omega_m0']
            self.params['Omega_l0'] = self.params['Omega_l0_ref']
            Omega_l0 = self.params['Omega_l0']

            phi_prime_arr = np.zeros(len(z_arr))
            E_arr = lcdm.compute_Ez_LCDM(z_arr, self.params['Omega_r0'], self.params['Omega_m0'])
            Omega_r_arr = lcdm.compute_Omega_r_z_LCDM(z_arr, self.params['Omega_r0'], self.params['Omega_m0'])
            Omega_m_arr = lcdm.compute_Omega_m_z_LCDM(z_arr, self.params['Omega_r0'], self.params['Omega_m0'])
            Omega_l_arr = lcdm.compute_Omega_l_z_LCDM(z_arr, self.params['Omega_r0'], self.params['Omega_m0'])
            
            E_prime_E_LCDM_arr = lcdm.compute_EprimeE_x_LCDM(x_arr, self.params['Omega_r0'], self.params['Omega_m0'])
            Omega_l_LCDM_arr = lcdm.compute_Omega_l_x_LCDM(x_arr, self.params['Omega_r0'], self.params['Omega_m0'])
            Omega_l_prime_LCDM_arr = lcdm.compute_Omega_l_prime_LCDM(E_prime_E_LCDM_arr, 1-self.params['Omega_r0']-self.params['Omega_m0'])

            E_prime_E_arr = np.copy(E_prime_E_LCDM_arr)
            E_prime_arr = E_prime_E_arr * E_arr

            phi_primeprime_arr = np.zeros(len(z_arr))

            A_arr = np.zeros(len(z_arr))

            Omega_phi_arr = np.zeros(len(z_arr))
            Omega_DE_arr = 1. - Omega_r_arr - Omega_m_arr

            Omega_phi_via_closure_arr = 1 - Omega_r_arr - Omega_m_arr - Omega_l_arr

            Omega_r_prime_arr = self.compute_Omega_r_prime(Omega_r_arr, E_arr, E_prime_arr)
            Omega_m_prime_arr = self.compute_Omega_m_prime(Omega_m_arr, E_arr, E_prime_arr)
            Omega_l_prime_arr = self.compute_Omega_l_prime(Omega_l_arr, E_arr, E_prime_arr)

            calB_arr = np.zeros(len(z_arr))
            calC_arr = np.zeros(len(z_arr))
            beta_arr = np.zeros(len(z_arr))

            chioverdelta_arr = np.zeros(len(z_arr))

            closure_variable = None
            roots_raw = None
            roots = None

            solver_success = None
        
        if forwards == False:

            # reverse direction of arrays...
            a_arr = a_arr[::-1]
            x_arr = x_arr[::-1]
            z_arr = z_arr[::-1]

            E_prime_E_LCDM_arr = E_prime_E_LCDM_arr[::-1]
            Omega_l_LCDM_arr = Omega_l_LCDM_arr[::-1]
            Omega_l_prime_LCDM_arr = Omega_l_prime_LCDM_arr[::-1]

            if roots is None and GR:
                E_arr = E_arr[::-1]
                E_prime_arr = E_prime_arr[::-1]
                E_prime_E_arr = E_prime_E_arr[::-1]
                phi_prime_arr = phi_prime_arr[::-1]
                phi_primeprime_arr = phi_prime_arr[::-1]
                Omega_m_arr = Omega_m_arr[::-1]
                Omega_r_arr = Omega_r_arr[::-1]
                Omega_l_arr = Omega_l_arr[::-1]
                Omega_phi_arr = Omega_phi_arr[::-1]
                Omega_phi_via_closure_arr = Omega_phi_via_closure_arr[::-1]
                Omega_DE_arr = Omega_DE_arr[::-1]
                Omega_m_prime_arr = Omega_m_prime_arr[::-1]
                Omega_r_prime_arr = Omega_r_prime_arr[::-1]
                Omega_l_prime_arr = Omega_l_prime_arr[::-1]
                A_arr = A_arr[::-1]
                calB_arr = calB_arr[::-1]
                calC_arr = calC_arr[::-1]
                beta_arr = beta_arr[::-1]
                chioverdelta_arr = chioverdelta_arr[::-1]

            else:
                E_arr = E_arr[:,::-1]
                E_prime_arr = E_prime_arr[:,::-1]
                E_prime_E_arr = E_prime_E_arr[:,::-1]
                phi_prime_arr = phi_prime_arr[:,::-1]
                phi_primeprime_arr = phi_prime_arr[:,::-1]
                Omega_m_arr = Omega_m_arr[:,::-1]
                Omega_r_arr = Omega_r_arr[:,::-1]
                Omega_l_arr = Omega_l_arr[:,::-1]
                Omega_phi_arr = Omega_phi_arr[:,::-1]
                Omega_phi_via_closure_arr = Omega_phi_via_closure_arr[:,::-1]
                Omega_DE_arr = Omega_DE_arr[:,::-1]
                Omega_m_prime_arr = Omega_m_prime_arr[:,::-1]
                Omega_r_prime_arr = Omega_r_prime_arr[:,::-1]
                Omega_l_prime_arr = Omega_l_prime_arr[:,::-1]
                A_arr = A_arr[:,::-1]
                calB_arr = calB_arr[:,::-1]
                calC_arr = calC_arr[:,::-1]
                beta_arr = beta_arr[:,::-1]
                chioverdelta_arr = chioverdelta_arr[:,::-1]
        
        E_like_arr = np.copy(E_arr)
        E_prime_like_arr = np.copy(E_prime_arr)

        if E_arr is not None:
            if np.isscalar(Omega_m0):
                E_arr /= E_like_arr[-1]
                E_prime_arr /= E_like_arr[-1]
            else:
                for idx in range(0, len(Omega_m0)):
                    E_arr[idx] /= E_like_arr[idx,-1]
                    E_prime_arr[idx] /= E_like_arr[idx,-1]

        self.output = {
            'a': a_arr,
            'x': x_arr,
            'z': z_arr,
            'E': E_arr,
            'E_like': E_like_arr,
            'H0': H0,
            'Omega_r0': Omega_r0,
            'Omega_m0': Omega_m0,
            'Omega_l0': Omega_l0,
            'E_prime': E_prime_arr,
            'E_prime_like': E_prime_like_arr,
            'E_prime/E': E_prime_E_arr,
            'E_prime/E_LCDM': E_prime_E_LCDM_arr,
            'phi_prime': phi_prime_arr,
            'phi_primeprime': phi_primeprime_arr,
            'Omega_m': Omega_m_arr,
            'Omega_r': Omega_r_arr,
            'Omega_l': Omega_l_arr,
            'Omega_l_LCDM': Omega_l_LCDM_arr,
            'Omega_phi': Omega_phi_arr,
            'Omega_phi_via_closure': Omega_phi_via_closure_arr,
            'Omega_DE': Omega_DE_arr,
            'Omega_m_prime': Omega_m_prime_arr,
            'Omega_r_prime': Omega_r_prime_arr,
            'Omega_l_prime': Omega_l_prime_arr,
            'Omega_l_prime_LCDM': Omega_l_prime_LCDM_arr,
            'A': A_arr,
            'calB': calB_arr,
            'calC': calC_arr,
            'beta': beta_arr,
            'chi/delta': chioverdelta_arr,
            'initialiser': {
                'z_start': z_start,
                'forwards': forwards,
                'closure_variable': closure_variable,
                'closure_values': roots,
                'closure_values_raw': roots_raw,
            },
            'solver_success': solver_success,
            'K_G3_G4_variables': [str(sym) for sym in self.sym['K_G3_G4_syms']],
            'K_G3_G4_values': self.params['K_G3_G4_values']
        }

        if compute_growth:
            self.get_linear_growth()
            self.get_linear_growth_2()

        return self.output

    
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
        
        if np.isscalar(self.output['Omega_m0']):

            fname_expansion = fname_prefix + 'expansion.txt'
            data = np.column_stack([self.output['a'], self.output['E'], self.output['E_prime/E']])
            np.savetxt(fname_expansion, data, fmt=['%.4e', '%.4e', '%.4e'])

            fname_force = fname_prefix + 'force.txt'
            data = np.column_stack([self.output['a'], self.output['chi/delta'], self.output['beta']])
            np.savetxt(fname_force, data, fmt=['%.4e', '%.4e', '%.4e'])
        
        else:
            
            for idx in range(0, len(self.output['Omega_m0'])):
                
                if len(self.output['Omega_m0']) != 1:
                    fname_expansion = fname_prefix + 'res_%i_expansion.txt' % idx
                else:
                    fname_expansion = fname_prefix + 'expansion.txt'
                data = np.column_stack([self.output['a'], self.output['E'][idx], self.output['E_prime/E'][idx]])
                np.savetxt(fname_expansion, data, fmt=['%.4e', '%.4e', '%.4e'])

                if len(self.output['Omega_m0']) != 1:
                    fname_force = fname_prefix + 'res_%i_force.txt' % idx
                else:
                    fname_force = fname_prefix + 'force.txt'
                data = np.column_stack([self.output['a'], self.output['chi/delta'][idx], self.output['beta'][idx]])
                np.savetxt(fname_force, data, fmt=['%.4e', '%.4e', '%.4e'])


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


    def clean(self):
        """
        Reinitialise the class. 
        """
        self.__init__()
