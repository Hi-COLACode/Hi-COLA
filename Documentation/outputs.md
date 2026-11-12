# Outputs and plotting

## Frontend tables

The Horndeski generator writes headerless ASCII tables (`%.4e`), starting at $a=1$ and running back to the earliest epoch:

| File | Columns |
| --- | --- |
| `*_expansion.txt` | $a$, $E=H/H_0$, $E'/E=d\ln H/d\ln a$ |
| `*_force.txt` | $a$, $\chi/\delta$, coupling $\beta$, $G_{G_4}/G_N$ |
| `*_full.txt` | $a$, $E$, $E'/E$, $\phi$, $\phi'$, $\Omega_m$, $\Omega_r$, $\Omega_\phi$, $\Omega_\Lambda$ (diagnostic only) |

The `HiCOLA` gravity model splines both tables against $a$ and uses $G_{\mathrm{eff}}/G_N=1+\beta$. The K-mouflage and DGP generators write three-column force tables, which the `kmouflage` gravity model reads. A quick structural check:

```python
import numpy as np
e = np.loadtxt("Hi-COLA_Output/cuGal_f0p2_14dp_WMAP_expansion.txt")
f = np.loadtxt("Hi-COLA_Output/cuGal_f0p2_14dp_WMAP_force.txt")
assert e.shape[1] == 3 and f.shape[1] == 4 and np.isfinite(e).all() and np.isfinite(f).all()
assert np.array_equal(e[:, 0], f[:, 0]) and (np.diff(e[:, 0]) < 0).all()
assert np.isclose(e[0, 0], 1) and np.isclose(e[0, 1], 1, rtol=1e-4)
assert e[-1, 0] <= 1 / (1 + 49)   # covers the backend's initial redshift
```

This does not show the solution is physically correct; compare with independent calculations and vary the resolution.

## Backend files

`output_folder`, `simulation_name`, and `output_redshifts` set where, what, and when. Snapshots use `output_particles` and `output_fileformat` (`GADGET` or `FML`); power spectra, halos, multipoles, and bispectra have their own toggles. Check each file's header.

**`cosmology_<name>.txt`** has nine columns, although its header names six: $a$, $H/H_0$, $d\ln H/d\ln a$, $\Omega_m$, photon radiation, relativistic neutrinos, massive-neutrino approximation, exact neutrino fraction, and a legacy dark-energy column. Columns 7 and 8 are alternative neutrino representations; do not add them. Column 9 is not the scalar's density; total dark energy is $1-\Omega_m-\Omega_r$.

**`gravitymodel_<name>_k<k>.txt`** starts with $a$, $G_{\mathrm{eff}}/G_N$, $D_1$, and a neutrino growth term at that $k$; later columns depend on the mode, so read [GravityModel::output](repo:HiCOLA/FML/COLASolver/src/GravityModel.h).

**`pofk*_z*.txt`** contains $k$ (h/Mpc), $P(k)$, and linear $P(k)$ ((Mpc/h)³). The header field `ShotnoiseSubtracted` follows `pofk_subtract_shotnoise`, which the template sets to true.

## Plotting tools

The packages install four commands, each taking the folder holding its input files:

| Command | Reads | Plots |
| --- | --- | --- |
| `hicola_plot_frontend` | `*_expansion.txt`, `*_force.txt` | $E$, $E'/E$, screening coefficient, coupling, $G_{G_4}/G_N$ |
| `hicola_plot_pofk` | `pofk*_z*.txt` | $P(k)$ at each output redshift |
| `hicola_plot_background` | `cosmology_*.txt` | $H/H_0$, $d\ln H/d\ln a$, density fractions |
| `hicola_plot_growth` | `gravitymodel_*_k*.txt` | $D_1$ and $G_{\mathrm{eff}}/G$ at each $k$ |

```sh
hicola_plot_pofk example-output/cubic_galileon/output --reference example-output/gr/output
```

`--reference` overlays a second run and adds ratio panels with bands set by `--bands` (default ±1% and ±2.5%; none if given no values). `--redshifts` and `--wavenumbers` select panels; `--label`, `--reference-label`, `--out`, and `--format` (png, pdf, svg) control the figure. Figures are saved to `plots/` beside an `output` folder; no window is opened. These figures are generated from the example simulations each time this documentation is built:

```{include} _generated/background.md
```

```{include} _generated/growth.md
```
