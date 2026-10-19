# Configure and run the backend

[Documentation home](index.md) · [Frontend](frontend.md) · [Outputs](outputs.md)

Build `nbody` using the [installation guide](installation.md) and generate both frontend tables first. Copy [Hi-COLA_params.lua](repo:HiCOLA/Parameter_files/Backend_parameter_files/Hi-COLA_params.lua) to your run directory, then edit that copy. The template contains placeholder paths and large production grids; it is not ready to execute unchanged.

## Connect the inputs

Set these Lua values, replacing paths with actual absolute paths:

```lua
cosmology_model = "HiCOLA"
gravity_model = "HiCOLA"
HiCOLA_expansion_filename = "/absolute/path/model_cosmology_expansion.txt"
HiCOLA_preforce_filename = "/absolute/path/model_cosmology_force.txt"
ic_input_filename = "/absolute/path/matterpower.dat"
output_folder = "/absolute/path/run-output"
```

The current [parameter reader](repo:HiCOLA/FML/COLASolver/src/ReadParameters.h) and [gravity class](repo:HiCOLA/FML/COLASolver/src/GravityModel_HiCOLA.h) require **`HiCOLA_preforce_filename`**. Replace the template's `HiCOLA_input_filename` assignment with that name. Table column one is scale factor `a`, despite template comments calling it `k`.

The initial spectrum is a separate input. For `ic_type_of_input = "powerspectrum"`, it contains `k` in `h/Mpc` and `P(k)` in `(Mpc/h)^3`; `ic_input_redshift` specifies its epoch. The bundled `LCDM_1306dot3219_matterpower_z0.000.dat` is an example, not a spectrum for every cosmology. Match spectrum, frontend densities, and backend `cosmology_OmegaCDM`, `cosmology_Omegab`, `cosmology_OmegaMNu`, and `cosmology_h`. The shipped frontend example is WMAP-based whereas the Lua template uses Planck values.

Inspect `ic_use_gravity_model_GR` and `ic_sigma8_normalization` when comparing modified-gravity and GR runs: they affect back-scaling and normalisation. Keep seeds and initial-condition choices recorded.

## Cosmology consistency and the dark-energy key

For the frontend's flat, massless-neutrino setup, use `cosmology_OmegaMNu = 0.0`. The backend's `cosmology_OmegaLambda` key represents **total present-day dark energy** in a Hi-COLA run, not just the true cosmological-constant component:

$$
\texttt{cosmology\_OmegaLambda}=\Omega_{\mathrm{DE},0}
=1-\Omega_{m0}-\Omega_{r0}.
$$

The frontend instead sets $\Omega_{\Lambda0}=(1-f_\phi)\Omega_{\mathrm{DE},0}$. Supplying that smaller cosmological-constant fraction to the backend would introduce an inconsistent curvature budget. The backend infers curvature from the density sum, whereas the frontend equations here assume spatial flatness.

Match the frontend's total radiation density to the backend radiation computed from `cosmology_TCMB_kelvin` and `cosmology_Neffective`, including relativistic neutrinos. Small rounding residuals may be treated as zero curvature by the backend; do not rely on that handling to conceal genuinely different cosmologies. See [initial conditions](initial-conditions.md) for spectrum backscaling and particle-snapshot input.

## Begin with small grids

For a workflow check, reduce the corresponding values in your Lua copy:

```lua
simulation_name = "hicola_smoke"
particle_Npart_1D = 16
ic_nmesh = particle_Npart_1D
force_nmesh = 32
pofk_nmesh = 32
timestep_nsteps = {5}
output_redshifts = {0.0}
output_particles = false
pofk = true
force_linear_massive_neutrinos = false
```

This checks execution and table wiring; these settings are not suitable for scientific results. Keep other required keys from the template because `all_parameters_must_be_in_file = true` makes even otherwise optional settings mandatory when read. Create the output directory before running.

Run from the repository root, substituting your Lua path:

```sh
OMP_NUM_THREADS=1 HiCOLA/FML/COLASolver/nbody /absolute/path/run.lua > run.log 2>&1
```

For an MPI build on a configured MPI system:

```sh
OMP_NUM_THREADS=1 mpiexec -n 2 HiCOLA/FML/COLASolver/nbody /absolute/path/run.lua > run.log 2>&1
```

Check the exit status immediately (`echo $?`), then inspect `run.log`, reported table paths, simulated epochs, and generated power-spectrum files. A successful small run should finish at the requested output redshift and produce finite output statistics. Investigate errors before increasing particle count, mesh sizes, MPI ranks, or threads. Grid memory grows roughly with the cube of mesh size, and multiple working grids add overhead.

## Screening

`gravity_model_screening = true` enables the Hi-COLA approximate screening prescription. Its settings include `gravity_model_HiCOLA_smoothing_filter`, `gravity_model_HiCOLA_smoothing_scale_over_boxsize`, `gravity_model_screening_enforce_largescale_linear`, and `gravity_model_screening_linear_scale_hmpc`. The smoothing scale is expressed relative to box size; the linear-transition scale is in `h/Mpc`. Disabling screening selects the linear fifth-force path.

Do not assume a new theory can use existing screening calibration unchanged. For K-mouflage applications, consult the model-specific paper linked in the [project README](repo:README.md); this page does not provide a validated K-mouflage preset. For general COLA settings and other gravity models, see the [bundled solver README](repo:HiCOLA/FML/COLASolver/README.md).
