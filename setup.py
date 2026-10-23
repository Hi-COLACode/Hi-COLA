from setuptools import setup, find_packages
setup(
    name = 'Hi_COLA',
    packages = find_packages(
    ),
    # Core runtime dependencies. These are imported by the HiCOLA Frontend
    # pipeline.
    # Ranges are chosen so CI resolves to the same major generation as a
    # modern local environment; relax or tighten together if that drifts.
    install_requires = [
        # Core numerical stack (array ops, ODE/quad solvers via scipy).
        # >=1.26 is the floor pandas 2.x/3.x and scipy 1.11+ build against.
        'numpy>=1.26',
        # Integrals/ODEs (scipy.integrate.solve_ivp, scipy.optimize.fsolve).
        # >=1.11 matches the numpy floor.
        'scipy>=1.11',
        # Symbolic algebra (sympy.lambdify) for building the model functions.
        'sympy>=1.12',
        # Plotting for diagnostics; >=3.8 is a stable modern floor.
        'matplotlib>=3.8',
        # INI-style config parsing (configobj.ConfigObj). 5.x is the current
        # stable line.
        'configobj>=5.0',
    ],
    extras_require = {
        # Test-only dependencies (installed via `pip install -e '.[tests]'`)
        'tests': [
            # Only imported in support.py for read_csv/concat
            'pandas>=2.0',
            # Test coverage reporting (coverage run -m unittest). >=7.6 is a
            # stable modern floor.
            'coverage>=7.6',
        ],
    },
)
