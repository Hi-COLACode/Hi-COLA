import numpy as np


def compute_Ez_LCDM(z, Omega_r0, Omega_m0):
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

    Returns
    -------
    E : float or array
        The normalised Hubble expansion rate.
    """
    # Note: replaces comp_E_LCDM
    Omega_L0 = 1. - Omega_m0 - Omega_r0
    E = np.sqrt(Omega_m0*(1+z)**3 + Omega_r0*(1+z)**4 + Omega_L0)
    return E


def compute_Hz_LCDM(z, Omega_r0, Omega_m0, H0):
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
    
    Returns
    -------
    H : float or array
        Hubble expansion rate.
    """
    # Note: replaces comp_H_LCDM
    H = H0 * compute_Ez_LCDM(z, Omega_r0, Omega_m0)
    return H


def compute_Omega_r_z_LCDM(z, Omega_r0, Omega_m0):
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
    
    Returns
    -------
    Omega_r : float or array
        Radiation fractional density.
    """
    # Note: replaces comp_Omega_r_LCDM
    E = compute_Ez_LCDM(z, Omega_r0, Omega_m0)
    Omega_r = Omega_r0 * (1+z)**4
    Omega_r /= E**2
    return Omega_r


def compute_Omega_m_z_LCDM(z, Omega_r0, Omega_m0):
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
    
    Returns
    -------
    Omega_m : float or array
        Matter fractional density.
    """
    # Note: replaces comp_Omega_m_LCDM
    E = compute_Ez_LCDM(z, Omega_r0, Omega_m0)
    Omega_m = Omega_m0 * (1+z)**3
    Omega_m /= E**2
    return Omega_m


def compute_Omega_l_z_LCDM(z, Omega_r0, Omega_m0):
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
    
    Returns
    -------
    Omega_L : float or array
        Lambda fractional density.
    """
    # Note: replaces comp_Omega_L_LCDM
    E = compute_Ez_LCDM(z, Omega_r0, Omega_m0)
    Omega_L0 = 1. - Omega_m0 - Omega_r0
    Omega_L = Omega_L0
    Omega_L /= E**2
    return Omega_L


def compute_EprimeE_x_LCDM(x, Omega_r0, Omega_m0):
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
    
    Returns
    -------
    EprimeE : float or array
        The derivative of the normalised Hubble function with respect to x, 
        divided by the normalised Hubble function.
    """
    # Note: replaces comp_E_prime_E_LCDM
    Omega_L0 = 1. - Omega_m0 - Omega_r0
    term1 = Omega_r0*np.exp(-4.*x) + Omega_m0*np.exp(-3.*x) + Omega_L0
    term2 = 4.*Omega_r0*np.exp(-4.*x) + 3.*Omega_m0*np.exp(-3.*x)
    EprimeE = -0.5*term2/term1
    return EprimeE


def compute_Omega_l_x_LCDM(x, Omega_r0, Omega_m0):
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
    
    Returns
    -------
    Omega_l : float or array
        The lambda (cosmological constant) fractional density.
    """
    Omega_l0 = 1. - Omega_m0 - Omega_r0
    term1 = Omega_r0*np.exp(-4.*x) + Omega_m0*np.exp(-3.*x) + Omega_l0
    Omega_l = Omega_l0/term1
    return Omega_l


def compute_Omega_l_prime_LCDM(EprimeE, Omega_l):
    """
    Compute the derivative of the fractional density for Lambda. Assuming a flat universe.

    Parameters
    ----------
    EprimeE : float or array
        The derivative of the normalised Hubble function with respect to x, 
        divided by the normalised Hubble function.
    Omega_l : float or array
        Lambda fractional density at equivalent times to EprimeE (Not to be confused for Omega_L0).
    
    Returns
    -------
    Omega_l_prime : float or array
    """
    # Note: replaces comp_Omega_DE_prime_LCDM
    Omega_l_prime = -2.*EprimeE*Omega_l
    return Omega_l_prime

