# Frontend configuration

[Documentation home](index.md) · [Quickstart](quickstart.md) · [Backend](backend.md)

The [generator](repo:HiCOLA/Frontend/generate_simulation_input.py) reads two INI files using [read_parameters.py](repo:HiCOLA/Frontend/read_parameters.py), builds symbolic field equations with [expression_builder.py](repo:HiCOLA/Frontend/expression_builder.py), and integrates them with [numerical_solver.py](repo:HiCOLA/Frontend/numerical_solver.py).

## Theory file

| Key | Purpose |
| --- | --- |
| `model_name` | Prefix of generated filenames. |
| `output_directory` | Destination, relative to the process working directory or absolute. |
| `K`, `G3`, `G4` | SymPy expressions for the reduced Horndeski functions. |
| Keys ending `_symbol` | Declare model constants, in parameter order. |
| Keys ending `_parameter` | Numerical values, in matching symbol order. |
| `f_phi` | Scalar fraction of today's dark-energy density; the reader sets `Omega_l0 = (1 - f_phi) * (1 - Omega_m0 - Omega_r0)`. |
| `Hubble0`, `phiprime0` | Present-day initial values in the model's normalisation. |
| `set_all_to_one` | Set all seven mass ratios to one. |
| `closure_guess_value` | Initial guess for the selected closure unknown. |

Write multiplication explicitly (`k_1*X`) and powers as `X**2`. `phi` and `X` are theory variables; other constants must be declared. A cubic Galileon definition is:

```ini
K = k_1*X
G3 = g_31*X
G4 = 0.5
k1_symbol = k_1
g31_symbol = g_31
k1_parameter = -1.4785008372906538
g31_parameter = 0.49283361243021795
```

This fragment is not a complete input file. Start with the [complete example](repo:HiCOLA/Parameter_files/Frontend_parameter_files/examples/cubic_Galileon/cuGal_f0p2.ini). Symbol and parameter lists are matched by file order, not by key name. Reordering only one list changes the model.

The supplied Galileon models use a rescaled Hubble function `U = E/E_dS`, with `E = H/H0`; their `Hubble0` need not equal one. The exporter divides the integrated Hubble function by the initial normalisation before writing `E`. Preserve the normalisation associated with your model parameters.

When `set_all_to_one = False`, the reader collects `M_` entries in file order. Use the template order: `M_pG4`, `M_KG4`, `M_G3s`, `M_sG4`, `M_G3G4`, `M_Ks`, `M_gp`. This less-tested path reads strings rather than explicitly converting all ratios to floats; validate it independently before relying on it.

## Closure equation

`use_constraint_eq_on_odeint_variables = True` selects an integration initial value; `False` selects a model parameter using zero-based `which_Horndeski_par`. The chosen value is replaced with `closure_guess_value` and solved using the closure equation.

For the current command-line reader, use `which_odeint_par = 1` for `phiprime0`. Index `0` selects the Hubble value but encounters a known undefined-variable defect in the closure solver. Although template comments list indices `2` and `3` for radiation and matter, the reader assigns into a two-element initial-condition list, so those indices fail before integration. The printed closure result identifies the actual solved value; the earlier printed parameter or initial value can be a guess.

## Numerical and cosmology file

| Key | Purpose |
| --- | --- |
| `cosmo_name` | Cosmology label in output filenames. |
| `Npoints` | Number of tabulated integration points. |
| `max_redshift` | Earliest redshift covered by the frontend. |
| `suppression_flag` | Suppress solver warnings; keep false while validating. |
| `GR_flag` | Select the frontend LCDM background path. |
| `threshold` | Threshold used by the numerical solver for its `A` coefficient. |
| `h` | Reduced Hubble constant. |
| `Omega_r0h2`, `Omega_b0h2`, `Omega_c0h2` | Physical radiation, baryon, and cold-dark-matter densities. |
| `Omega_r0`, `Omega_m0` | Direct density inputs, or literal `None` to derive from physical densities. |

With `Omega_m0 = None`, matter is `(Omega_b0h2 + Omega_c0h2)/h**2`. With `Omega_r0 = None`, radiation is `Omega_r0h2/h**2`. **Keep `Omega_m0 = None` in the current reader:** its direct-matter branch mistakenly reads `Omega_r0`. Setting a direct matter value therefore does not supply that value correctly.

Choose `max_redshift` to cover the backend's initial redshift and all output epochs. Validate convergence by varying resolution and numerical settings; successful integration alone does not establish physical viability or screening accuracy. For theory conventions and model derivations, consult the [equation and normalisation reference](equations.md).
