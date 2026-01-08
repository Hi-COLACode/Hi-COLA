import numpy as np
import time

from . import redshift, utils


class StandardModel:

    """
    A class for computing the LCDM background expansion and growth functions to construct 
    non-linear cosmological simulations using the Hi-COLA Backend.
    """
    
    # Initialising the main class

    def __init__(self):
        """
        Initialises the standard LCDM model class.
        """

        # Enable print statements
        self.verbose = True

        # Parameter values
        self.params = {}
        
        # Output dictionary
        self.outputs = {}

        # Timer dictionary
        self.timer = {'t0': None}

        # Constants dictionary
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

        # Interpolation functions
        self.interp = {}
    

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
    
    # Densities for components of the Universe.

    ## Photons are given by the CMB temperature

    def _get_C_gamma(self):
        """
        Retrieve C_gamma constant.
        """
        self.const['C_gamma'] =  8 * (np.pi**3) * self.const['G[eV]'] * (self.params['Tcmb0'] * self.const['kB[eV]'])**4
        self.const['C_gamma'] /= 45*(self.const['H0unit[eV]']**2)
    
    
    def get_rho_g(self, a, h):
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
    

    def get_rho_g_prime(self, a, h):
        """
        Returns the E^2 Omega_gamma (photon) density derivative.

        Parameters
        ----------
        a : float or array
            Scale factor.
        h : float
            Little h = H0 * 1e-2.
        """
        return -4*self.get_rho_g(a, h)
    

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
        return self.get_rho_g(a, h)/(E**2)

    
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


    def get_rho_nu(self, a, h, mnu):
        """
        Returns the neutrino density, i.e. Omega_nu*E**2.

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
        rho_nu : float or array
            Neutrino density, i.e. Omega_nu*E**2.
        """
        if self.neutrino_table is None:
            self._tabulate_IJy()
        y = a * mnu / (self.params['Tnu0']*self.const['kB[eV]'])
        if utils.isscalar(y):
            if y <= self.neutrino_table['ymin']:
                rho_nu = self.neutrino_table['Iy_interp_low'](y)
            elif y >= self.neutrino_table['ymax']:
                rho_nu = self.neutrino_table['Iy_interp_high'](y)
            else:
                rho_nu = self.neutrino_table['Iy_interp'](y)
        else:
            rho_nu = np.zeros(len(y))
            cond = np.where((y <= self.neutrino_table['ymin']))[0]
            rho_nu[cond] = self.neutrino_table['Iy_interp_low'](y[cond])
            cond = np.where((y >= self.neutrino_table['ymax']))[0]
            rho_nu[cond] = self.neutrino_table['Iy_interp_high'](y[cond])
            cond = np.where((y > self.neutrino_table['ymin']) & (y < self.neutrino_table['ymax']))[0]
            rho_nu[cond] = self.neutrino_table['Iy_interp'](y[cond])
        rho_nu *= self.const['C_nu']/((h**2) * (a**4))
        return rho_nu
    
    
    def get_rho_nu_prime(self, a, h, mnu):
        """
        Returns the neutrino density (i.e. Omega_nu*E**2) prime.

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
        rho_nu_prime : float or array
            Neutrino density derivative wrt log a.
        """
        rho_nu_prime = -4*self.get_rho_nu(a, h, mnu)
        y = a * mnu / (self.params['Tnu0']*self.const['kB[eV]'])
        if utils.isscalar(y):
            if y <= self.neutrino_table['ymin']:
                _rho_nu_prime = self.neutrino_table['Iy_prime_interp_low'](y)
            elif y >= self.neutrino_table['ymax']:
                _rho_nu_prime = self.neutrino_table['Iy_prime_interp_high'](y)
            else:
                _rho_nu_prime = self.neutrino_table['Iy_prime_interp'](y)
        else:
            _rho_nu_prime = np.zeros(len(y))
            cond = np.where((y <= self.neutrino_table['ymin']))[0]
            _rho_nu_prime[cond] = self.neutrino_table['Iy_prime_interp_low'](y[cond])
            cond = np.where((y >= self.neutrino_table['ymax']))[0]
            _rho_nu_prime[cond] = self.neutrino_table['Iy_prime_interp_high'](y[cond])
            cond = np.where((y > self.neutrino_table['ymin']) & (y < self.neutrino_table['ymax']))[0]
            _rho_nu_prime[cond] = self.neutrino_table['Iy_prime_interp'](y[cond])
        _rho_nu_prime *= self.const['C_nu']/((h**2) * (a**4))
        rho_nu_prime += _rho_nu_prime
        return rho_nu_prime
    
    
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
        return self.get_rho_nu(a, h, mnu)/(E**2)


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
        if utils.isscalar(y):
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
        if utils.isscalar(y):
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
    
    def get_rho_b(self, a, Omega_b0):
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
    

    def get_rho_b_prime(self, a, Omega_b0):
        """
        Returns the baryon density derivative.

        Parameters
        ----------
        a : float or array
            Scale factor.
        Omega_b0 : float
            Baryon fractional density at redshift zero.
        """
        return -3*self.get_rho_b(a, Omega_b0)
    

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
        return self.get_rho_b(a, Omega_b0)/(E**2)
    

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


    def get_rho_c(self, a, Omega_c0):
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
    

    def get_rho_c_prime(self, a, Omega_c0):
        """
        Returns the cold dark matter density derivative.

        Parameters
        ----------
        a : float or array
            Scale factor.
        Omega_c0 : float
            Cold dark matter fractional density at redshift zero.
        """
        return -3*self.get_rho_c(a, Omega_c0)
    

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
        return self.get_rho_c(a, Omega_c0)/(E**2)
    

    def compute_Omega_c_prime(self, Omega_c, E, E_prime):
        """
        Computes Omega cold dark matter prime.

        Parameters
        ----------
        Omega_c : float or array
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

    def get_rho_l(self, a, Omega_l0, w0=-1., wa=0.):
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
    

    def get_rho_l_prime(self, a, Omega_l0, w0=-1., wa=0.):
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
            return -3*(1+w0+wa*(1-a))*self.get_rho_l(a, Omega_l0, w0=w0, wa=wa)
    

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
        return self.get_rho_l(a, Omega_l0, w0=w0, wa=wa)/(E**2)
    

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

        E2 = self.get_rho_g(a, h)
        E2 += self.get_rho_c(a, Omega_c0)
        E2 += self.get_rho_b(a, Omega_b0)
        E2 += self.get_rho_l(a, Omega_l0, w0=w0, wa=wa)
        E2 += self.get_rho_nu(a, h, mnu[0]) + self.get_rho_nu(a, h, mnu[1]) + self.get_rho_nu(a, h, mnu[2])

        # Sqrt to obtain E

        E = np.sqrt(E2)

        return E
    

    def compute_E_prime_LCDM(self, a, h, Omega_b0, Omega_c0, mnu=[0.,0.,0.], w0=-1., wa=0.):
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
        E_prime : float or array
            The derivative of normalised Hubble expansion rate.
        """

        E = self.compute_E_LCDM(a, h, Omega_b0, Omega_c0, mnu=mnu, w0=-1., wa=0.)

        # Obtain photon and neutrino fractional density today to obtain dark energy density
        # assuming a spatially flat universe.

        Omega_g0 = self.get_Omega_g(1., h, 1.)
        Omega_nu0 = self.get_Omega_nu(1., h, mnu[0], 1.) + self.get_Omega_nu(1., h, mnu[1], 1.) + self.get_Omega_nu(1., h, mnu[2], 1.)
        
        Omega_l0 = 1. - Omega_g0 - Omega_nu0 - Omega_c0 - Omega_b0

        E_prime = self.get_rho_g_prime(a, h)
        E_prime += self.get_rho_nu_prime(a, h, mnu[0]) + self.get_rho_nu_prime(a, h, mnu[1]) + self.get_rho_nu_prime(a, h, mnu[2])
        E_prime += self.get_rho_b_prime(a, Omega_b0)
        E_prime += self.get_rho_c_prime(a, Omega_c0)
        E_prime += self.get_rho_l_prime(a, Omega_l0, w0=w0, wa=wa)

        E_prime /= 2*E

        return E_prime


    def set_cosmo_params(self, H0_ref, Omega_c0_ref, Omega_b0_ref, w0=-1., wa=0., Tcmb=2.7255, Tnu0=1.9518, mnu=[0.,0.,0.], Neff=3.04):
        """
        Set cosmological parameters.

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
        self.params['Neff'] = Neff
        self.params['Omega_nu10_ref'] = self.get_Omega_nu(1., 1e-2*self.params['H0_ref'], self.params['mnu'][0], 1.) 
        self.params['Omega_nu20_ref'] = self.get_Omega_nu(1., 1e-2*self.params['H0_ref'], self.params['mnu'][1], 1.) 
        self.params['Omega_nu30_ref'] = self.get_Omega_nu(1., 1e-2*self.params['H0_ref'], self.params['mnu'][2], 1.)
        self.params['Omega_l0_LCDM'] = 1. - self.params['Omega_g0_ref'] - self.params['Omega_nu10_ref'] - self.params['Omega_nu20_ref'] - self.params['Omega_nu30_ref']- self.params['Omega_c0_ref'] - self.params['Omega_b0_ref']
        self.params['Omega_l0_ref'] = self.params['Omega_l0_LCDM']
        self.params['w0'] = w0
        self.params['wa'] = wa
    

    def comp_w_eff(self):
        """
        Computes effective equation of state.
        """
        if self.output['E'] is None or self.output['E_prime'] is None:
            self.output['w_eff'] = None
        else:
            self.output['w_eff'] = - 1 - (2/3)*self.output['E_prime']/self.output['E']
    

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
                    Sigma = np.copy(mu)
                    gamma = np.ones(len(mu))
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
                        Sigma[idx] = np.copy(mu[idx])
                        gamma[idx] = np.ones(len(mu[idx]))
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
    

    def get_r_drag(self):
        """
        Compute r_drag given in Mpc using the approximation equation 2 of 2503.14738.
        """

        if self.output['Omega_b0'] is None or self.output['Omega_c0'] is None or self.output['H0'] is None:
            
            r_drag = None
        
        else:

            if utils.isscalar(self.output['Omega_b0']):
                
                w_b = self.output['Omega_b0']*(1e-2*self.output['H0'])**2
                w_bc = (self.output['Omega_b0'] + self.output['Omega_c0'])*(1e-2*self.output['H0'])**2

                r_drag = 147.05
                r_drag *= (w_b/0.02236)**-0.13
                r_drag *= (w_bc/0.1432)**-0.23
                r_drag *= (self.params['Neff']/3.04)**-0.1
                
            else:

                r_drag = np.zeros(np.shape(self.output['Omega_b0']))

                for idx in range(0, len(self.output['Omega_b0'])):

                    w_b = self.output['Omega_b0'][idx]*(1e-2*self.output['H0'][idx])**2
                    w_bc = (self.output['Omega_b0'][idx] + self.output['Omega_c0'][idx])*(1e-2*self.output['H0'][idx])**2

                    r_drag[idx] = 147.05
                    r_drag[idx] *= (w_b/0.02236)**-0.13
                    r_drag[idx] *= (w_bc/0.1432)**-0.23
                    r_drag[idx] *= (self.params['Neff']/3.04)**-0.1
        
        self.output['r_drag'] = r_drag
    

    def get_a_eq(self):
        """
        Computes the radiation-matter equality scale factor, using the Hu-Sugiyama fitting function.
        """
        if self.output['Omega_b0'] is None or self.output['Omega_c0'] is None or self.output['H0'] is None:
            
            a_eq = None
        
        else:

            if utils.isscalar(self.output['Omega_b0']):
                
                fnu = 0.405
                Theta27 = 2.7255/2.7
                w_bc = (self.output['Omega_b0'] + self.output['Omega_c0'])*(1e-2*self.output['H0'])**2
                a_eq = 2.35*1e-5*Theta27**4
                a_eq /= w_bc*(1-fnu)

                ## Note: below is the previous equation which computes this more exactly, however if we are
                ## going to use the Hu-Sugiyama approximation we need to do so self-consistently.
                # a_eq = self.const['C_gamma'] + self.const['C_nu']*7.*(np.pi**4)/40.
                # a_eq /= (self.output['Omega_b0'] + self.output['Omega_c0'])*(1e-2*self.output['H0'])**2
                
            else:

                a_eq = np.zeros(np.shape(self.output['Omega_b0']))

                for idx in range(0, len(self.output['Omega_b0'])):

                    fnu = 0.405
                    Theta27 = 2.7255/2.7
                    w_bc = (self.output['Omega_b0'][idx] + self.output['Omega_c0'][idx])*(1e-2*self.output['H0'][idx])**2
                    a_eq[idx] = 2.35*1e-5*Theta27**4
                    a_eq[idx] /= w_bc*(1-fnu)

                    # a_eq[idx] = self.const['C_gamma'] + self.const['C_nu']*7.*(np.pi**4)/40.
                    # a_eq[idx] /= (self.output['Omega_b0'][idx] + self.output['Omega_c0'][idx])*(1e-2*self.output['H0'][idx])**2
        
        self.output['a_eq'] = a_eq
    

    def get_z_star(self):
        """
        Computes the redshift for last scattering, using the Hu-Sugiyama fitting function.
        """
        if self.output['Omega_b0'] is None or self.output['Omega_c0'] is None or self.output['H0'] is None:
            
            z_star = None
        
        else:

            if utils.isscalar(self.output['Omega_b0']):
                
                w_b = self.output['Omega_b0']*(1e-2*self.output['H0'])**2
                w_bc = (self.output['Omega_b0']+self.output['Omega_c0'])*(1e-2*self.output['H0'])**2

                g1 = 0.0783*(w_b**-0.238)/(1 + 39.5*w_b**0.763)
                g2 = 0.560/(1+21.1*w_b**1.81)

                z_star = 1048 * (1 + 0.00124*w_b**-0.738)*(1 + g1*w_bc**g2)
                
            else:

                z_star = np.zeros(np.shape(self.output['Omega_b0']))

                for idx in range(0, len(self.output['Omega_b0'])):

                    w_b = self.output['Omega_b0'][idx]*(1e-2*self.output['H0'][idx])**2
                    w_bc = (self.output['Omega_b0'][idx]+self.output['Omega_c0'][idx])*(1e-2*self.output['H0'][idx])**2

                    g1 = 0.0783*(w_b**-0.238)/(1 + 39.5*w_b**0.763)
                    g2 = 0.560/(1+21.1*w_b**1.81)

                    z_star[idx] = 1048. * (1 + 0.00124*w_b**-0.738)*(1 + g1*w_bc**g2)

        self.output['a_star'] = redshift.z2a(z_star)
    

    def get_r_star(self):
        """
        Compute r_star.
        """
        if self.output['Omega_b0'] is None or self.output['Omega_c0'] is None or self.output['H0'] is None or self.output['a_eq'] is None or self.output['a_star'] is None:
            
            r_star = None
        
        else:

            if utils.isscalar(self.output['Omega_b0']):
                
                w_b = self.output['Omega_b0']*(1e-2*self.output['H0'])**2
                Omega_bc = self.output['Omega_b0']+self.output['Omega_c0']
                Theta27 = 2.726/2.7

                H0 = self.output['H0'] # in km s^-1 Mpc^-1
                H0 /= self.const['c[km/s]'] # in Mpc^-1

                ## Note this is the exact equation but we need to use the Hu-Sugiyama approximations
                ## self-consistently
                # R = (3/4)*self.get_rho_b(self.output['a_star'], self.output['Omega_b0'])
                # R /= self.get_rho_g(self.output['a_star'], self.output['H0']*1e-2)

                # R_eq = (3/4)*self.get_rho_b(self.output['a_eq'], self.output['Omega_b0'])
                # R_eq /= self.get_rho_g(self.output['a_eq'], self.output['H0']*1e-2)

                R = 31.5*w_b*(Theta27**-4)*1e3/redshift.a2z(self.output['a_star'])
                R_eq = 31.5*w_b*(Theta27**-4)*1e3/redshift.a2z(self.output['a_eq'])
                
                r_star = 2*np.sqrt(3)/3
                r_star /= np.sqrt(Omega_bc*H0**2)
                r_star *= np.sqrt(self.output['a_eq']/R_eq)
                r_star *= np.log((np.sqrt(1+R) + np.sqrt(R+R_eq))/(1+np.sqrt(R_eq)))
                
            else:

                r_star = np.zeros(np.shape(self.output['Omega_b0']))

                for idx in range(0, len(self.output['Omega_b0'])):

                    Omega_bc = self.output['Omega_b0'][idx]+self.output['Omega_c0'][idx]

                    w_b = self.output['Omega_b0'][idx]*(1e-2*self.output['H0'][idx])**2
                    Omega_bc = self.output['Omega_b0'][idx]+self.output['Omega_c0'][idx]
                    Theta27 = 2.726/2.7

                    H0 = self.output['H0'][idx] # in km s^-1 Mpc^-1
                    H0 /= self.const['c[km/s]'] # in Mpc^-1
                    
                    # R = (3/4)*self.get_rho_b(self.output['a_star'][idx], self.output['Omega_b0'][idx])
                    # R /= self.get_rho_g(self.output['a_star'][idx], self.output['H0'][idx]*1e-2)

                    # R_eq = (3/4)*self.get_rho_b(self.output['a_eq'][idx], self.output['Omega_b0'][idx])
                    # R_eq /= self.get_rho_g(self.output['a_eq'][idx], self.output['H0'][idx]*1e-2)

                    R = 31.5*w_b*(Theta27**-4)*1e3/redshift.a2z(self.output['a_star'][idx])
                    R_eq = 31.5*w_b*(Theta27**-4)*1e3/redshift.a2z(self.output['a_eq'][idx])
                    
                    r_star[idx] = 2*np.sqrt(3)/3
                    r_star[idx] /= np.sqrt(Omega_bc*H0**2)
                    r_star[idx] *= np.sqrt(self.output['a_eq'][idx]/R_eq)
                    r_star[idx] *= np.log((np.sqrt(1+R) + np.sqrt(R+R_eq))/(1+np.sqrt(R_eq)))

        self.output['r_star'] = r_star


    def get_theta_star(self):
        """
        Computes the theta_star and produces comoving distance interpolating function.
        """
        if self.output['Omega_b0'] is None or self.output['Omega_c0'] is None or self.output['H0'] is None or self.output['a_eq'] is None or self.output['a_star'] is None or self.output['r_star'] is None:
            
            theta_star = None
        
        else:
            
            from scipy.interpolate import interp1d

            if utils.isscalar(self.output['Omega_b0']):
                
                self.interp['Dc%i_vs_x[Mpc]' % 0] = interp1d(self.output['x'], self.output['Dc']/(self.output['H0']*1e-2), kind='linear')

                theta_star = self.output['r_star']/self.interp['Dc%i_vs_x[Mpc]' % 0](redshift.a2x(self.output['a_star']))
                
            else:

                theta_star = np.zeros(np.shape(self.output['Omega_b0']))

                for idx in range(0, len(self.output['Omega_b0'])):

                    self.interp['Dc%i_vs_x[Mpc]' % idx] = interp1d(self.output['x'], self.output['Dc'][idx]/(self.output['H0'][idx]*1e-2), kind='linear')

                    theta_star[idx] = self.output['r_star'][idx]/self.interp['Dc%i_vs_x[Mpc]' % idx](redshift.a2x(self.output['a_star'][idx]))
            
        self.output['theta_star'] = theta_star

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

        z_start = z_arr[0]

        self.output['initialiser'] = {}
        self.output['initialiser']['z_start'] = z_start
        self.output['initialiser']['forwards'] = forwards


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

        self.output['Ehat'] = None
        self.output['Ehat_prime'] = None
        self.output['phihat'] = None
        self.output['phihat_prime'] = None
        self.output['phihat_primeprime'] = None

        self.output['E'] = self.compute_E_LCDM(self.output['a'], 1e-2*self.params['H0'], self.params['Omega_b0'], self.params['Omega_c0'], self.params['mnu'], w0=self.params['w0'], wa=self.params['wa'])
        self.output['E_prime'] = self.compute_E_prime_LCDM(self.output['a'], 1e-2*self.params['H0'], self.params['Omega_b0'], self.params['Omega_c0'], self.params['mnu'], w0=self.params['w0'], wa=self.params['wa'],)
        self.output['phi'] = np.zeros(len(self.output['z']))
        self.output['phi_prime'] = np.zeros(len(self.output['z']))
        self.output['phi_primeprime'] = np.zeros(len(self.output['z']))

        self.output['rhohat_phi'] = None
        self.output['rhohat_g'] = None
        self.output['rhohat_b'] = None
        self.output['rhohat_c'] = None
        self.output['rhohat_l'] = None
        self.output['rhohat_nu1'] = None
        self.output['rhohat_nu2'] = None
        self.output['rhohat_nu3'] = None

        self.output['rho_phi'] = np.zeros(len(self.output['z']))
        self.output['rho_g'] = self.get_rho_g(self.output['a'], 1e-2*self.params['H0'])
        self.output['rho_b'] = self.get_rho_b(self.output['a'], self.params['Omega_b0'])
        self.output['rho_c'] = self.get_rho_c(self.output['a'], self.params['Omega_c0'])
        self.output['rho_l'] = self.get_rho_l(self.output['a'], self.params['Omega_l0'], w0=self.params['w0'], wa=self.params['wa'])
        self.output['rho_nu1'] = self.get_rho_nu(self.output['a'], 1e-2*self.params['H0'], self.params['mnu'][0])
        self.output['rho_nu2'] = self.get_rho_nu(self.output['a'], 1e-2*self.params['H0'], self.params['mnu'][1])
        self.output['rho_nu3'] = self.get_rho_nu(self.output['a'], 1e-2*self.params['H0'], self.params['mnu'][2])

        self.output['Omega_phi'] = np.zeros(len(self.output['z']))
        self.output['Omega_g'] = self.get_Omega_g(self.output['a'], 1e-2*self.params['H0'], self.output['E'])
        self.output['Omega_b'] = self.get_Omega_b(self.output['a'], self.params['Omega_b0'], self.output['E'])
        self.output['Omega_c'] = self.get_Omega_c(self.output['a'], self.params['Omega_c0'], self.output['E'])
        self.output['Omega_l'] = self.get_Omega_l(self.output['a'], self.params['Omega_l0'], self.output['E'], w0=self.params['w0'], wa=self.params['wa'])
        self.output['Omega_nu1'] = self.get_Omega_nu(self.output['a'], 1e-2*self.params['H0'], self.params['mnu'][0], self.output['E'])
        self.output['Omega_nu2'] = self.get_Omega_nu(self.output['a'], 1e-2*self.params['H0'], self.params['mnu'][1], self.output['E'])
        self.output['Omega_nu3'] = self.get_Omega_nu(self.output['a'], 1e-2*self.params['H0'], self.params['mnu'][2], self.output['E'])

        self.output['w_l'] = self.compute_w_l(self.output['a'], w0=self.params['w0'], wa=self.params['wa'])
        self.output['w_nu1'] = self.compute_w_nu(self.output['a'], self.params['mnu'][0])
        self.output['w_nu2'] = self.compute_w_nu(self.output['a'], self.params['mnu'][1])
        self.output['w_nu3'] = self.compute_w_nu(self.output['a'], self.params['mnu'][2])


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
        w_l_arr = self.output['w_l']

        G_G_4_G_N = np.ones(len(x_arr))

        E_prime_arr = self.compute_E_prime_LCDM(a_arr, self.params['H0']*1e-2, self.params['Omega_b0'], self.params['Omega_c0'], mnu=self.params['mnu'], w0=self.params['w0'], wa=self.params['wa'])

        Omega_l_prime_arr = self.compute_Omega_l_prime(a_arr, Omega_l_arr, E_arr, E_prime_arr)

        A_arr = np.zeros(len(x_arr))

        Omega_phi_arr = self.output['Omega_phi']
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

        M_star_sq_arr = np.ones(len(x_arr))
        alpha_M_arr = np.zeros(len(x_arr))
        alpha_B_arr = np.zeros(len(x_arr))
        alpha_B_prime_arr = np.zeros(len(x_arr))
        alpha_K_arr = np.zeros(len(x_arr))
        
        tilde_calE_arr = np.zeros(len(x_arr))
        tilde_calP_arr = np.zeros(len(x_arr))
        w_phi_arr = np.zeros(len(x_arr))

        w_DE_arr = w_l_arr*Omega_l_arr + w_phi_arr*Omega_phi_arr
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
        self.output['A'] = A_arr
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
        self.output['tilde_calE'] = tilde_calE_arr
        self.output['tilde_calP'] = tilde_calP_arr
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
                'a', 'z', 'x', 
                'Ehat', 'Ehat_prime', 'phihat', 'phihat_prime', 'phihat_primeprime', 
                'E', 'E_prime', 'phi', 'phi_prime', 'phi_primeprime', 
                'rhohat_phi', 'rhohat_g', 'rhohat_b', 'rhohat_c', 'rhohat_l', 'rhohat_nu1', 'rhohat_nu2', 'rhohat_nu3', 
                'rho_phi', 'rho_g', 'rho_b', 'rho_c', 'rho_l', 'rho_nu1', 'rho_nu2', 'rho_nu3', 
                'Omega_phi', 'Omega_g', 'Omega_b', 'Omega_c', 'Omega_l', 'Omega_nu1', 'Omega_nu2', 'Omega_nu3', 
                'w_nu1', 'w_nu2', 'w_nu3', 'w_l',
                # derived quantities
                'H', 'Dc', 'G_G_4/G_N',
                'A', 'Omega_DE', 'w_DE', 'Omega_phi_via_closure', 
                'Omega_g_prime', 'Omega_b_prime', 'Omega_c_prime', 'Omega_l_prime', 
                'Omega_nu1_prime', 'Omega_nu2_prime', 'Omega_nu3_prime', 
                'calB', 'calC', 'beta', 'chi/delta', 'M_star_sq',
                'alpha_M', 'alpha_B', 'alpha_B_prime', 'alpha_K' ,
                'tilde_calE', 'tilde_calP', 'w_phi', 'D', 'Q_s', 'c_s_sq_D', 'c_s_sq'
            ]
        else:
            keys = [
                'a', 'z', 'x', 
                'Ehat', 'Ehat_prime', 'phihat', 'phihat_prime', 'phihat_primeprime', 
                'E', 'E_prime', 'phi', 'phi_prime', 'phi_primeprime', 
                'rhohat_g', 'rhohat_b', 'rhohat_c', 'rhohat_l', 'rhohat_nu1', 'rhohat_nu2', 'rhohat_nu3', 
                'rho_g', 'rho_b', 'rho_c', 'rho_l', 'rho_nu1', 'rho_nu2', 'rho_nu3', 
                'Omega_g', 'Omega_b', 'Omega_c', 'Omega_l', 'Omega_nu1', 'Omega_nu2', 'Omega_nu3', 
                'w_nu1', 'w_nu2', 'w_nu3', 'w_l',
            ]

        for key in keys:
            if self.output[key] is not None:
                if self.output[key].ndim == 1:
                    self.output[key] = self.output[key][::-1]
                elif self.output[key].ndim == 2:
                    self.output[key] = self.output[key][:,::-1]


    def run_solver(
            self, z_max=2000., Npoints=1000, forwards=True, derived=True
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
        derived : bool, optional
            Instructs the solver whether derived quantities should be computed.
        
        Returns
        -------
        outputs : dict
            Dictionary containing numerical solver solutions.
        """

        if self.verbose:
            print('Hi-COLA: computing LCDM expansion and linear perturbation')

        self.output = {}

        self._run_solver_start(z_max, Npoints, forwards)
        
        if self.verbose:
            print(' - Running in GR mode!')

        if self.verbose:
            print(' - Computing expansion history.')
        
        self._run_solver_ODE_GR()

        if derived:
            if self.verbose:
                print(' - Computing main derived quantities.')
            
            self._run_solver_derived_GR()
        
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
                data = np.column_stack([self.output['a'], self.output['E'][idx], self.output['E_prime'][idx]/self.output['E'][idx]])
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
