# Testing

[Documentation home](index.md) · [Installation](installation.md)

## Unit tests

From a source checkout:

```sh
python -m pip install -e '.[tests]'
MPLBACKEND=Agg python -m unittest discover -s tests -v
```

Measure Python statement and branch coverage with:

```sh
MPLBACKEND=Agg python -m coverage run -m unittest discover -s tests -v
python -m coverage report
python -m coverage html
```

`htmlcov/index.html` shows uncovered lines and branches. The suite tests the Python frontend, plotting tools, and the Python backend launcher; it does not measure C++ execution or collect anything under `HiCOLA/FML`. Continuous integration runs it on Python 3.11 to 3.15 and uploads coverage reports. [tests/README.md](repo:tests/README.md) describes the reference values, tolerances, tests marked as expected failures for known defects, and Codecov setup.

## Regression tests

The **Regression tests** workflow runs the [small example simulations](examples.md) for GR, cubic Galileon, and K-mouflage on a pinned baseline revision and on the candidate revision, then compares them:

- **Regression summary statistics** compares matter power spectra, expansion and growth tables, particle counts, mean velocities, and RMS speeds against explicit tolerances. Make it a required check.
- **Regression particle reproducibility** compares particle positions, velocities, and force grids. It is advisory.

The `regression-summary-plots` artifact contains the [`hicola_plot_*`](plotting.md) figures for every case of both revisions, `overview.png` showing every group of checks against its tolerance, and velocity statistics for each case. Particle-level tolerances are defined in `tests/regression_tests/compare.py`.

## Package builds

The **Build** workflow builds the pip wheel and source archive and the native conda packages, installs them, and checks them on Linux, macOS, and Windows. See [RELEASING.md](repo:RELEASING.md) for publication.
