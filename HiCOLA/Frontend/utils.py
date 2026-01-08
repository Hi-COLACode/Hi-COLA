import numpy as np


def isscalar(x):
    """More general isscalar function to prevent 0 dimensional numpy arrays
    from being misidentified as arrays even though they are actually scalar
    variables.
    """
    if type(x).__module__ == np.__name__:
        if len(x.shape) == 0:
            return True
        else:
            return False
    else:
        return np.isscalar(x)
