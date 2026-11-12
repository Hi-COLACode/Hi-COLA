# Frontend configuration

The [generator](repo:HiCOLA/Frontend/generate_simulation_input.py) reads a theory file and a numerical/cosmology file, builds the field equations symbolically, solves the closure condition, and integrates the background. Start from the [complete example](repo:HiCOLA/Parameter_files/Frontend_parameter_files/examples/cubic_Galileon/cuGal_f0p2.ini); the blank `horndeski_parameters.ini` template cannot run until `K`, `G3`, and `G4` are filled in.

## Theory file

| Key | Purpose |
| --- | --- |
| `model_name`, `output_directory` | Output filename prefix and destination (relative to the working directory). |
| `K`, `G3`, `G4` | SymPy expressions in `phi`, `X`, and declared constants. |
| `*_symbol`, `*_parameter` | Constants and their values, **matched by file order**, not by name. |
| `f_phi` | Scalar fraction of present dark energy: $\Omega_{\Lambda0}=(1-f_\phi)(1-\Omega_{m0}-\Omega_{r0})$. |
| `Hubble0`, `phiprime0` | Present-day initial values in the model's normalisation. |
| `set_all_to_one` | Set all seven [mass ratios](theory.md#dimensions-and-mass-ratios) to one. |
| `closure_guess_value` | Starting guess for the closure unknown. |

```ini
K = k_1*X
G3 = g_31*X
G4 = 0.5
k1_symbol = k_1
g31_symbol = g_31
k1_parameter = -1.4785008372906538
g31_parameter = 0.49283361243021795
```

Write products and powers explicitly (`k_1*X`, `X**2`). The supplied Galileon models use the [rescaled Hubble function](theory.md#initial-conditions-and-rescaling) $U=E/E_{\mathrm{dS}}$, so `Hubble0` need not be one; the exporter divides by it before writing $E$. With `set_all_to_one = False`, `M_` entries are read in the template's order (`M_pG4`, `M_KG4`, `M_G3s`, `M_sG4`, `M_G3G4`, `M_Ks`, `M_gp`) as strings; validate this less-tested path before relying on it.

## Closure

`use_constraint_eq_on_odeint_variables = True` solves for an initial value selected by `which_odeint_par`; `False` solves for model parameter `which_Horndeski_par` (zero-based). Use `which_odeint_par = 1` (`phiprime0`): index `0` fails on an undefined variable, and indices `2` and `3` fail because only two initial values are stored. The printed closure result is the solved value; a value printed earlier may be the guess.

## Numerical and cosmology file

| Key | Purpose |
| --- | --- |
| `cosmo_name` | Cosmology label in output filenames. |
| `Npoints`, `max_redshift` | Number of output points and earliest redshift; cover the backend's initial redshift. |
| `suppression_flag` | Suppress solver warnings; keep `False` while validating. |
| `GR_flag` | Use the ΛCDM background path. |
| `threshold` | Threshold for the solver's $A$ coefficient. |
| `h`, `Omega_r0h2`, `Omega_b0h2`, `Omega_c0h2` | Reduced Hubble constant and physical densities. |
| `Omega_r0`, `Omega_m0` | Direct densities, or `None` to derive them from the physical densities. |

**Keep `Omega_m0 = None`:** the direct-matter branch reads `Omega_r0` by mistake. Successful integration does not establish physical viability; check convergence by varying `Npoints` and the numerical settings.
