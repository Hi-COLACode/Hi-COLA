import numpy as np


def compute_Ez_LCDM(z, Omega_r0, Omega_m0, w0=-1., wa=0., mnu=None):
    """
    Compute the normalised Hubble function E=H/H0 for LCDM. Assuming a flat universe.

    Parameters
    ----------
    z : float or array
        Redshift.
    Omega_r0 : float
        Radiation fractional density today (i.e. z=0).
    Omega_m0 : float
        Matter fractional density today (i.e. z=0).
    w0 : float, optional
        Set dark energy equation of state today.
    wa : float, optional
        Set dark energy equation of state gradient.
    mnu : list, optional
        Sets the neutrino mass.

    Returns
    -------
    E : float or array
        The normalised Hubble expansion rate.
    """
    # Note: replaces comp_E_LCDM
    Omega_L0 = 1. - Omega_m0 - Omega_r0
    if w0 == -1. and wa == 0.:
        E = np.sqrt(Omega_m0*(1+z)**3 + Omega_r0*(1+z)**4 + Omega_L0)
    else:
        a = 1/(1+z)
        E = np.sqrt(Omega_m0*(1+z)**3 + Omega_r0*(1+z)**4 + Omega_L0*(a**(-3*(1+w0+wa)))*np.exp(3*wa*(a-1)))
    return E


def compute_Hz_LCDM(z, Omega_r0, Omega_m0, H0, w0=-1., wa=0.):
    """
    Compute the Hubble function for LCDM. Assuming a flat universe.

    Parameters
    ----------
    z : float or array
        Redshift.
    Omega_r0 : float
        Radiation fractional density today (i.e. z=0).
    Omega_m0 : float
        Matter fractional density today (i.e. z=0).
    H0 : float
        Hubble constant
    w0 : float, optional
        Set dark energy equation of state today.
    wa : float, optional
        Set dark energy equation of state gradient.
    
    Returns
    -------
    H : float or array
        Hubble expansion rate.
    """
    # Note: replaces comp_H_LCDM
    H = H0 * compute_Ez_LCDM(z, Omega_r0, Omega_m0, w0=w0, wa=wa)
    return H


def compute_Omega_r_z_LCDM(z, Omega_r0, Omega_m0, w0=-1., wa=0.):
    """
    Compute the radiation fractional density at a given redshift in LCDM. Assuming a flat universe.

    Parameters
    ----------
    z : float or array
        Redshift.
    Omega_r0 : float
        Radiation fractional density today (i.e. z=0).
    Omega_m0 : float
        Matter fractional density today (i.e. z=0).
    w0 : float, optional
        Set dark energy equation of state today.
    wa : float, optional
        Set dark energy equation of state gradient.
    
    Returns
    -------
    Omega_r : float or array
        Radiation fractional density.
    """
    # Note: replaces comp_Omega_r_LCDM
    E = compute_Ez_LCDM(z, Omega_r0, Omega_m0, w0=w0, wa=wa)
    Omega_r = Omega_r0 * (1+z)**4
    Omega_r /= E**2
    return Omega_r


def compute_Omega_m_z_LCDM(z, Omega_r0, Omega_m0, w0=-1., wa=0.):
    """
    Compute the matter fractional density at a given redshift in LCDM. Assuming a flat universe.

    Parameters
    ----------
    z : float or array
        Redshift.
    Omega_r0 : float
        Radiation fractional density today (i.e. z=0).
    Omega_m0 : float
        Matter fractional density today (i.e. z=0).
    w0 : float, optional
        Set dark energy equation of state today.
    wa : float, optional
        Set dark energy equation of state gradient.
    
    Returns
    -------
    Omega_m : float or array
        Matter fractional density.
    """
    # Note: replaces comp_Omega_m_LCDM
    E = compute_Ez_LCDM(z, Omega_r0, Omega_m0, w0=w0, wa=wa)
    Omega_m = Omega_m0 * (1+z)**3
    Omega_m /= E**2
    return Omega_m


def compute_Omega_l_z_LCDM(z, Omega_r0, Omega_m0, w0=-1., wa=0.):
    """
    Compute the fractional density for Lambda at a given redshift in LCDM. Assuming a flat universe.

    Parameters
    ----------
    z : float or array
        Redshift.
    Omega_r0 : float
        Radiation fractional density today (i.e. z=0).
    Omega_m0 : float
        Matter fractional density today (i.e. z=0).
    w0 : float, optional
        Set dark energy equation of state today.
    wa : float, optional
        Set dark energy equation of state gradient.
    
    Returns
    -------
    Omega_L : float or array
        Lambda fractional density.
    """
    # Note: replaces comp_Omega_L_LCDM
    E = compute_Ez_LCDM(z, Omega_r0, Omega_m0, w0=w0, wa=wa)
    Omega_L0 = 1. - Omega_m0 - Omega_r0
    if w0 == -1. and wa == 0.:
        Omega_L = Omega_L0
    else:
        a = 1/(1+z)
        Omega_L = Omega_L0*(a**(-3*(1+w0+wa)))*np.exp(3*wa*(a-1))
    Omega_L /= E**2
    return Omega_L


def compute_EprimeE_x_LCDM(x, Omega_r0, Omega_m0, w0=-1., wa=0.):
    """
    Computes Eprime/E as a function of x = log(a).

    Parameters
    ----------
    x : float or array
        The natural logarithm of the scale factor a.
    Omega_r0 : float
        Radiation fractional density today (i.e. z=0).
    Omega_m0 : float
        Matter fractional density today (i.e. z=0).
    w0 : float, optional
        Set dark energy equation of state today.
    wa : float, optional
        Set dark energy equation of state gradient.
    
    Returns
    -------
    EprimeE : float or array
        The derivative of the normalised Hubble function with respect to x, 
        divided by the normalised Hubble function.
    """
    # Note: replaces comp_E_prime_E_LCDM
    Omega_L0 = 1. - Omega_m0 - Omega_r0
    if w0 == -1. and wa == 0.:
        term1 = Omega_r0*np.exp(-4.*x) + Omega_m0*np.exp(-3.*x) + Omega_L0
        term2 = 4.*Omega_r0*np.exp(-4.*x) + 3.*Omega_m0*np.exp(-3.*x)
        EprimeE = -0.5*term2/term1
    else:
        a = np.exp(x)
        term1 = Omega_r0*np.exp(-4.*x) + Omega_m0*np.exp(-3.*x) + Omega_L0*(a**(-3*(1+w0+wa)))*np.exp(3*wa*(a-1))
        term2 = 4.*Omega_r0*np.exp(-4.*x) + 3.*Omega_m0*np.exp(-3.*x) + (3*(1+w0+wa) + 3*wa*a)*Omega_L0*(a**(-3*(1+w0+wa)))*np.exp(3*wa*(a-1))
        EprimeE = -0.5*term2/term1
    return EprimeE


def compute_Omega_l_x_LCDM(x, Omega_r0, Omega_m0, w0=-1., wa=0.):
    # TODO: Remove, just convert x to a -> a=np.exp(x)
    """
    Computes Eprime/E as a function of x = log(a).

    Parameters
    ----------
    x : float or array
        The natural logarithm of the scale factor a.
    Omega_r0 : float
        Radiation fractional density today (i.e. z=0).
    Omega_m0 : float
        Matter fractional density today (i.e. z=0).
    w0 : float, optional
        Set dark energy equation of state today.
    wa : float, optional
        Set dark energy equation of state gradient.
    
    Returns
    -------
    Omega_l : float or array
        The lambda (cosmological constant) fractional density.
    """
    Omega_l0 = 1. - Omega_m0 - Omega_r0
    if w0 == -1. and wa == 0.:
        term1 = Omega_r0*np.exp(-4.*x) + Omega_m0*np.exp(-3.*x) + Omega_l0
        Omega_l = Omega_l0/term1
    else:
        a = np.exp(x)
        term1 = Omega_r0*np.exp(-4.*x) + Omega_m0*np.exp(-3.*x) + Omega_l0*(a**(-3*(1+w0+wa)))*np.exp(3*wa*(a-1))
        Omega_l = Omega_l0*(a**(-3*(1+w0+wa)))*np.exp(3*wa*(a-1))/term1
    return Omega_l


def compute_Omega_l_prime_LCDM(x, EprimeE, Omega_l, w0=-1., wa=0.):
    """
    Compute the derivative of the fractional density for Lambda. Assuming a flat universe.

    Parameters
    ----------
    x : float or array
        The natural logarithm of the scale factor a.
    EprimeE : float or array
        The derivative of the normalised Hubble function with respect to x, 
        divided by the normalised Hubble function.
    Omega_l : float or array
        Lambda fractional density at equivalent times to EprimeE (Not to be confused for Omega_L0).
    w0 : float, optional
        Set dark energy equation of state today.
    wa : float, optional
        Set dark energy equation of state gradient.
    
    Returns
    -------
    Omega_l_prime : float or array
    """
    # Note: replaces comp_Omega_DE_prime_LCDM
    if w0 == -1. and wa == 0.:
        Omega_l_prime = -2.*EprimeE*Omega_l
    else:
        a = np.exp(x)
        Omega_l_prime = (-3*(1+w0+wa) + 3*wa*a - 2*EprimeE)*Omega_l
    return Omega_l_prime

