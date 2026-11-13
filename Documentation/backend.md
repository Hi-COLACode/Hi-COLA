# Backend configuration

Copy the [Lua template](repo:HiCOLA/Parameter_files/Backend_parameter_files/Hi-COLA_params.lua) and edit the copy; it has placeholder paths and production-sized grids. The [example simulations](quickstart.md#run-the-example-simulations) are smaller, complete alternatives. `all_parameters_must_be_in_file = true` makes every key read for the selected modes mandatory, so keep the template's keys.

## Inputs

```lua
cosmology_model = "HiCOLA"
gravity_model = "HiCOLA"
HiCOLA_expansion_filename = "/path/model_cosmology_expansion.txt"
HiCOLA_preforce_filename = "/path/model_cosmology_force.txt"
ic_input_filename = "/path/matterpower.dat"
output_folder = "/path/run-output"
```

The tables' first column is scale factor, although template comments call it `k`. Match the backend cosmology (`cosmology_OmegaCDM`, `cosmology_Omegab`, `cosmology_h`, and radiation from `cosmology_TCMB_kelvin` and `cosmology_Neffective`) to the frontend; the shipped frontend example is WMAP-based and the Lua template Planck-based. Use `cosmology_OmegaMNu = 0.0`.

`cosmology_OmegaLambda` is the **total** present dark energy, $1-\Omega_{m0}-\Omega_{r0}$, not the frontend's $\Omega_{\Lambda0}=(1-f_\phi)\Omega_{\mathrm{DE},0}$. The backend infers curvature from the density sum, so supplying the smaller value would add spurious curvature.

## Initial conditions

With `ic_type_of_input = "powerspectrum"`, the input file has one header line, then $k$ in h/Mpc and $P(k)$ in (Mpc/h)³ at redshift `ic_input_redshift`. It is rescaled to `ic_initial_redshift` with the backend growth factor,

$$
P_{\mathrm{ini}}(k)=\left[\frac{D_1(k,a_{\mathrm{ini}})}{D_1(k,a_{\mathrm{input}})}\right]^2P_{\mathrm{input}}(k),
$$

using GR growth if `ic_use_gravity_model_GR = true`. `ic_sigma8_normalization` can also renormalise it. Ideally use a linear Einstein–Boltzmann spectrum for the same theory. Otherwise, an approximation is $P^{\mathrm{MG}}(k,z_{\mathrm{ini}})=[D_1^{\mathrm{MG}}(k,z_{\mathrm{ini}})/D_1^{\Lambda\mathrm{CDM}}(k,z_{\mathrm{target}})]^2P^{\Lambda\mathrm{CDM}}(k,z_{\mathrm{target}})$, with compatible early-time normalisations. Set `ic_input_redshift` to the file's actual epoch so a spectrum is not rescaled twice. Transfer-function inputs (`transferfunction`, `transferinfofile`) are shown in the [solver examples](repo:HiCOLA/FML/COLASolver/examples/example_run_different_input).

To start from a GADGET snapshot, set `ic_random_field_type = "read_particles"` and `ic_reconstruct_gadgetfilepath` to the file prefix (usually ending `/gadget`). The LPT fields are reconstructed from the particles using `ic_reconstruct_assigment_method` (spelled so in the code), `ic_reconstruct_interlacing`, and `ic_reconstruct_smoothing_filter`; match `ic_nmesh` to the particle grid.

Gaussian realisations differ only by `ic_random_seed` (record `ic_random_generator` too). For **fixed-and-paired** runs, keep the seed, set `ic_fix_amplitude = true`, and run twice with `ic_reverse_phases = false` and `true`; average the two spectra. This reduces variance in the mean but does not give the Gaussian covariance.

## Resolution

| Setting | Production example | Note |
| --- | --- | --- |
| `simulation_boxsize` | 512 or 1024 | Box side in Mpc/h. |
| `particle_Npart_1D` | 512 or 1024 | Per axis; the total is its cube. |
| `force_nmesh` | `3 * particle_Npart_1D` | Three cells per interparticle spacing. |
| `timestep_nsteps` | `{50}` | With `timestep_scalefactor_spacing = "linear"`. |

These follow the cited [COLA accuracy study](references.md#cola-growth-and-simulation-choices) and are starting points, not optimal values. Particle mass is

$$
m_{\mathrm{part}}\simeq8.325\times10^{10}
\left(\frac{\Omega_{m0}}{0.3}\right)\left(\frac{L}{1000\,h^{-1}\mathrm{Mpc}}\right)^3
\left(\frac{1000}{N_{\mathrm{1D}}}\right)^3h^{-1}M_\odot .
$$

Small scales are limited by the number of time steps, the force cell $L/N_{\mathrm{force}}$, particle shot noise, and the model-dependent screening and LPT approximations; the box limits large scales. Establish the usable $k$ range by convergence tests. For GR/modified-gravity ratios, keep seeds, normalisation, box, meshes, and output redshifts identical. Memory grows with the cube of the mesh and particle counts.

## Screening

`gravity_model_screening = true` enables the approximate screening of [the theory page](theory.md#force-and-growth); `false` uses the linear fifth force. It is tuned by `gravity_model_HiCOLA_smoothing_filter`, `gravity_model_HiCOLA_smoothing_scale_over_boxsize`, `gravity_model_screening_enforce_largescale_linear`, and `gravity_model_screening_linear_scale_hmpc` (h/Mpc). Validate the calibration for each new theory; K-mouflage needs the settings of its [paper](references.md#hi-cola-applications-and-validation). Other settings are described in the [solver README](repo:HiCOLA/FML/COLASolver/README.md).

## Run the simulation

For a first check, reduce `particle_Npart_1D` to 16, the meshes to 32, and `timestep_nsteps` to `{5}`, and set `output_particles = false`. Create the output folder, then run:

```sh
OMP_NUM_THREADS=1 mpiexec -n 2 HiCOLA/FML/COLASolver/nbody run.lua > run.log 2>&1
```

With the conda package or a native build, run from Python instead:

```python
from HiCOLA import backend
backend.run("run.lua", cwd="/path/to/simulation")
```

Relative Lua paths resolve from `cwd`. A failed run raises `subprocess.CalledProcessError`; `capture_output` and `timeout` are passed to `subprocess.run`. `backend.run` starts a single process; to use an MPI Makefile build, set `HICOLA_BACKEND_EXECUTABLE` and launch with `mpiexec` yourself. EuclidEmulator2 predictions need its optional FML header, which is not bundled.

Check the exit status and `run.log`, then the [outputs](outputs.md). Save particles (`output_particles = true`, `output_redshifts`, `output_fileformat`) only when you need them later, for example to reuse exact initial conditions or build an external halo catalogue.
