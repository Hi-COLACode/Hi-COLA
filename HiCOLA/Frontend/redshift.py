import numpy as np


def z2a(z):
    """
    Converts redshift to scale factor a.

    Parameters
    ----------
    z : float or array_like
        Redshift.
    """
    return 1/(1+z)


def a2z(a):
    """
    Converts scale factor to redshift.

    Parameters
    ----------
    a : float or array_like
        scale factor
    """
    return 1./a  - 1.


def a2x(a):
    """
    Converts scale factor a to x=log(a)

    Parameters
    ----------
    a : float or array_like
        scale factor
    """
    return np.log(a)


def x2a(x):
    """
    Converts the log of the scale factor to scale factor.

    Parameters
    ----------
    x : float or array_like
        log of the scale factor.
    """
    return np.exp(x)


def z2x(z):
    """
    Converts redshift to the log of the scale factor a.

    Parameters
    ----------
    z : float or array_like
        Redshift.
    """
    return a2x(z2a(z))


def x2z(x):
    """
    Converts the log of the scale factor to redshift.

    Parameters
    ----------
    x : float or array_like
        log of the scale factor.
    """
    return a2z(x2a(x))