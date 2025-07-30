import numpy as np
import sympy as sym


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


    def get_Pde(self):
        """
        UNSURE what this is computing, Pressure for dark energy?
        """
        # Note: replaces Pde function
        self.sym['H'] = sym.Symbol('H')
        self.sym['Hprime'] = sym.Symbol('Hprime')
        self.sym['Meff'] = sym.Symbol('M_eff')
        self.sym['Mp'] = sym.Symbol('Mp')
        self.sym['Ms'] = sym.Symbol('Ms')
        # The following phi's are defined as Tilde phi, however they cal G3 and G4 which rely on phi, so this must just be normal phi, called something else for some reason?
        # We will assume it's the same thing and not define it again here as a separate variable.
        # self.sym['Tildephi'] = sym.Symbol('Tildephi')
        # self.sym['Tildephiprime'] = sym.Symbol('Tildephiprime')
        # self.sym['Tildephiprimeprime'] = sym.Symbol('Tildephiprimeprime')
        Mfrac = (self.sym['Mp']**2)/(self.sym['Meff']**2)
        term1 = 2*self.sym['X']*(self.symfunc['G3phi'] + self.sym['H']*(self.sym['Hprime']*self.sym['Ms']*self.sym['phiprime']
            + self.sym['H']*self.sym['Ms']*self.sym['phiprimeprime'])*self.symfunc['G3x'])
        term2 = 2*self.symfunc['G4phi']*(self.sym['H']*(self.sym['Hprime']*self.sym['Ms']*self.sym['phiprime'] 
            + self.sym['H']*self.sym['Ms']*self.sym['phiprimeprime']) + 2*(self.sym['H']**2)*self.sym['Ms']*self.sym['phiprime'])
        term3 = 4*self.sym['X']*self.symfunc['G4phiphi']
        return Mfrac*(self.symfunc['K'] - term1 + term2 + term3) + ((self.sym['H']**2)*self.sym['Omega_r']*(self.sym['Mp']**2))* (Mfrac - 1)


    def get_rhode(self):
        """
        UNSURE if this is a defunct function...
        """
        # the original function takes in phidot and phidotdot and appeared to be missing a * sign. I think this is a dead or unused function.
        # Note: replaces rhode function
        # might need to reinitialise the H, Hprime, etc from get_Pde, we will assume they are already defined.
        Mfrac = (self.sym['Mp']**2)/(self.sym['Meff']**2)
        term1 = 2*self.sym['X']*self.symfunc['Kx'] - self.symfunc['K'] + 6*self.sym['X']*self.sym['Ms']*self.sym['phiprime']*(self.sym['H']**2)*self.symfunc['G3x']
        term2 = 2*self.sym['X']*self.symfunc['G3phi'] + 6*(self.sym['H']**2)*self.sym['Ms']*self.sym['phiprime']*self.symfunc['G4phi']
        return (3*(self.sym['H']**2)*(self.sym['Mp']**2))*(self.sym['Omega_r'] + self.sym['Omega_m'])*(Mfrac - 1) + Mfrac*(term1 - term2)
    

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


    def get_coupling_factor(self):
        """
        The coupling in equation 3.13 in https://arxiv.org/abs/2209.01666, i.e. the deviation from 1.
        """
        self.get_alpha1()
        self.get_alpha2()
        self.get_calC()
        self.symfunc['coupling'] = -1.*(self.symfunc['alpha1'] + self.symfunc['alpha2'])*self.symfunc['calC']


    def construct_model(self):
        """
        Constructs Horndeski model with user defined functions.
        """

        if self._check_symfunc_keys(['K', 'G3', 'G4']) == False:
            assert False, 'Functions for K, G3 and G4 remain undefined.'
        else:
            self.get_K_derivatives()
            self.get_G3_derivatives()
            self.get_G4_derivatives()

            Xreal = 0.5*(self.sym['E']**2.)*self.sym['phiprime']**2.

            self.get_EprimeE()
            EprimeE = self.symfunc['EprimeE'].subs(self.sym['X'], Xreal)

            self.get_EprimeE_safe()
            EprimeE_safe = self.symfunc['EprimeE_safe'].subs(self.sym['X'], Xreal)

            self.get_phiprimeprime()
            phiprimeprime = self.symfunc['phiprimeprime'].subs(self.sym['X'], Xreal)

            self.get_phiprimeprime_safe()
            phiprimeprime_safe = self.symfunc['phiprimeprime_safe'].subs(self.sym['X'], Xreal)

            self.get_A()
            A_func = self.symfunc['A'].subs(self.sym['X'], Xreal)

            self.get_B2()
            B2_func = self.symfunc['B2'].subs(self.sym['X'], Xreal)

            self.get_Omega_phi()
            Omega_phi = self.symfunc['Omega_phi'].subs(self.sym['X'], Xreal)

            self.get_fried_closure()
            fried_closure = self.symfunc['fried_closure'].subs(self.sym['X'], Xreal)

            self.get_alpha0()
            alpha0_func = self.symfunc['alpha0'].subs(self.sym['X'], Xreal)

            self.get_alpha1()
            alpha1_func = self.symfunc['alpha1'].subs(self.sym['X'], Xreal)

            self.get_alpha2()
            alpha2_func = self.symfunc['alpha2'].subs(self.sym['X'], Xreal)

            self.get_beta0()
            beta0_func = self.symfunc['beta0'].subs(self.sym['X'], Xreal)

            self.get_calB()
            calB_func = self.symfunc['calB'].subs(self.sym['X'], Xreal)

            self.get_calC()
            calC_func = self.symfunc['calC'].subs(self.sym['X'], Xreal)

            self.get_coupling_factor()
            coupling_fac = self.symfunc['coupling'].subs(self.sym['X'], Xreal)

            # Lambdify functions
            self.lambda_funcs['B2_lambda'] = sym.lambdify([self.sym['E'], self.sym['phiprime'], *self.sym['K_G3_G4_syms']], B2_func, "scipy")
            self.lambda_funcs['fried_closure_lambda'] = sym.lambdify(
                [self.sym['E'], self.sym['phiprime'], self.sym['Omega_r'], self.sym['Omega_m'], self.sym['Omega_l'], *self.sym['K_G3_G4_syms']],
                fried_closure
            )
            self.lambda_funcs['EprimeE_lambda'] = sym.lambdify([self.sym['E'], self.sym['phiprime'], self.sym['Omega_r'], self.sym['Omega_l'], *self.sym['K_G3_G4_syms']], EprimeE, "scipy")
            self.lambda_funcs['EprimeE_safe_lambda'] = sym.lambdify(
                [self.sym['E'], self.sym['phiprime'], self.sym['Omega_r'], self.sym['Omega_l'], self.sym['threshold'], self.sym['threshold_sign'], *self.sym['K_G3_G4_syms']],
                EprimeE_safe, "scipy"
            )
            self.lambda_funcs['phiprimeprime_lambda'] = sym.lambdify([self.sym['E'], self.sym['Eprime'], self.sym['phiprime'], *self.sym['K_G3_G4_syms']], phiprimeprime, "scipy")
            self.lambda_funcs['phiprimeprime_safe_lambda'] = sym.lambdify(
                [self.sym['E'], self.sym['Eprime'], self.sym['phiprime'], self.sym['threshold'], self.sym['threshold_sign'], *self.sym['K_G3_G4_syms']], 
                phiprimeprime_safe, "scipy"
            )
            self.lambda_funcs['Omega_phi_lambda'] = sym.lambdify([self.sym['E'], self.sym['phiprime'], *self.sym['K_G3_G4_syms']], Omega_phi)
            self.lambda_funcs['A_lambda'] = sym.lambdify([self.sym['E'], self.sym['phiprime'], *self.sym['K_G3_G4_syms']], A_func, "scipy")
            self.lambda_funcs['alpha0_lambda'] = sym.lambdify([self.sym['E'], self.sym['Eprime'], self.sym['phiprime'], self.sym['phiprimeprime'], *self.sym['K_G3_G4_syms']], alpha0_func)
            self.lambda_funcs['alpha1_lambda'] = sym.lambdify([self.sym['E'], self.sym['Eprime'], self.sym['phiprime'], self.sym['phiprimeprime'], *self.sym['K_G3_G4_syms']], alpha1_func)
            self.lambda_funcs['alpha2_lambda'] = sym.lambdify([self.sym['E'], self.sym['Eprime'], self.sym['phiprime'], self.sym['phiprimeprime'], *self.sym['K_G3_G4_syms']], alpha2_func)
            self.lambda_funcs['beta0_lambda'] = sym.lambdify([self.sym['E'], self.sym['Eprime'], self.sym['phiprime'], self.sym['phiprimeprime'], *self.sym['K_G3_G4_syms']], beta0_func)
            self.lambda_funcs['calB_lambda'] = sym.lambdify([self.sym['E'], self.sym['Eprime'], self.sym['phiprime'], self.sym['phiprimeprime'], *self.sym['K_G3_G4_syms']], calB_func)
            self.lambda_funcs['calC_lambda'] = sym.lambdify([self.sym['E'], self.sym['Eprime'], self.sym['phiprime'], self.sym['phiprimeprime'], *self.sym['K_G3_G4_syms']], calC_func)
            self.lambda_funcs['coupling_fac_lambda'] = sym.lambdify([self.sym['E'], self.sym['Eprime'], self.sym['phiprime'], self.sym['phiprimeprime'], *self.sym['K_G3_G4_syms']], coupling_fac)

    
    def run_solver(self):
        pass


    def clean(self):
        self.__init__()
