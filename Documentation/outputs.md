# Outputs and validation

[Documentation home](index.md) · [Backend](backend.md) · [Troubleshooting](troubleshooting.md)

## Frontend tables

The generator writes whitespace-separated ASCII files without headers. Both have three columns and use `%.4e` formatting (five significant digits). The writer reverses the integration arrays to store increasing scale factor.

| File | Column 1 | Column 2 | Column 3 |
| --- | --- | --- | --- |
| `*_expansion.txt` | `a` | `E = H/H0` | `d ln H / d ln a` |
| `*_force.txt` | `a` | `chi / delta` | Coupling factor |

The backend splines each table against `a`. Its linear effective Newton constant is `Geff/G = 1 + coupling`; `chi/delta` controls the approximate screening function. These tables are not functions of Fourier wavenumber.

Filenames include `model_name` and `cosmo_name`. If either target exists, both names are passed through the renaming helper; inspect the printed names and update the Lua input paths. A fresh output directory per run avoids ambiguous pairs and the helper's suffix limitations.

Run this structural check in your Python environment, replacing the two paths and initial redshift:

```python
import numpy as np

expansion = np.loadtxt("Hi-COLA_Output/cuGal_f0p2_14dp_WMAP_expansion.txt")
force = np.loadtxt("Hi-COLA_Output/cuGal_f0p2_14dp_WMAP_force.txt")
initial_redshift = 49.0
for table in (expansion, force):
    assert table.ndim == 2 and table.shape[1] == 3
    assert np.isfinite(table).all()
    assert (np.diff(table[:, 0]) > 0).all()
    assert table[0, 0] > 0
    assert table[0, 0] <= 1 / (1 + initial_redshift)
    assert np.isclose(table[-1, 0], 1.0)
assert np.array_equal(expansion[:, 0], force[:, 0])
assert (expansion[:, 1] > 0).all()
assert np.isclose(expansion[-1, 1], 1.0, rtol=1e-4)
```

These checks catch malformed tables and missing time coverage. They do not establish that the equations, chosen solution branch, or screening approximation are physically correct. Compare background histories and growth to independent references, and repeat with numerical resolution changes.

## Backend products

`output_folder` selects the destination and `simulation_name` labels the run. `output_redshifts` selects output epochs. Particle snapshots are controlled by `output_particles` and `output_fileformat` (`GADGET` or `FML`). Power spectra use the `pofk` settings; halo finding, multipoles, and bispectra have separate toggles in the Lua template.

Check each output's own headers and the [solver documentation](repo:HiCOLA/FML/COLASolver/README.md) before interpreting columns. Preserve the simulation log alongside the output. An exit status of zero is necessary, but also verify the requested epochs, expected files, finite statistics, and intended cosmology and gravity model.

For comparisons, match initial spectrum, random seed, phase/amplitude settings, box size, meshes, and time stepping. Establish convergence and statistical tolerances before treating differences in matter power spectra as physical effects.

## Cosmology diagnostic table

`cosmology_<simulation_name>.txt` is written by [Cosmology::output](repo:HiCOLA/FML/COLASolver/src/Cosmology.h). Its current rows contain **nine columns**, even though the header lists only six:

| Column | Quantity |
| --- | --- |
| 1 | Scale factor $a$ |
| 2 | $H/H_0$ |
| 3 | $d\ln H/d\ln a$ |
| 4 | $\Omega_m(a)$ |
| 5 | Photon radiation fraction returned by `get_OmegaR` |
| 6 | Relativistic-neutrino fraction returned by `get_OmegaNu` |
| 7 | Massive-neutrino matter approximation returned by `get_OmegaMNu` |
| 8 | Neutrino energy fraction returned by `get_OmegaNu_exact` |
| 9 | Legacy dark-energy quantity returned by `get_OmegaLambda` |

The tabulated expansion is splined from the frontend table; density diagnostics depend on consistent backend inputs. In a Hi-COLA run, column 9 is not generally the true cosmological-constant fraction and does not by itself describe an evolving scalar. For the supported massless-neutrino, flat setup, total dark energy can be reconstructed as $1-\Omega_m-\Omega_r$, with total radiation formed from photon and relativistic-neutrino contributions. Do not add the approximate and exact neutrino diagnostic columns together: they are alternative representations, not independent species.

## Gravity and growth diagnostics

`gravitymodel_<simulation_name>_k<value>.txt` files describe gravity and growth at the indicated fixed Fourier scale. The first columns are $a$, $G_{\mathrm{eff}}/G_N$, $D_1$, and the neutrino growth diagnostic. Transfer-data growth columns may then be inserted when transfer data are present, followed by second/third-order growth and source terms. Inspect the header and [GravityModel::output](repo:HiCOLA/FML/COLASolver/src/GravityModel.h) rather than assuming the remaining columns have a single layout for every mode. The first three columns are usually the starting point for comparing expansion-driven growth and coupling.

## Power-spectrum interpretation

The ordinary matter-spectrum output identifies $k$ in `h/Mpc`, nonlinear $P(k)$, and a linear reference in `(Mpc/h)^3`. Massive-neutrino modes can write additional component spectra. The header records `ShotnoiseSubtracted`; check that field and the explicit `pofk_subtract_shotnoise` setting rather than assuming all spectra include or exclude shot noise. The supplied Lua template sets subtraction true, while other defaults or examples may differ.

Use the [resolution guide](resolution.md) to decide which scales are trustworthy and estimate particle masses. Use the [paired-fixed recipe](initial-conditions.md#fixed-and-paired-runs) when comparing mean spectra with reduced sampling variance. The choice to save snapshots depends on whether later particle-level analysis or exact initial-condition reuse is required.
