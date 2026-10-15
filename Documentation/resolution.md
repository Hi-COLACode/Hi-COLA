# Accuracy, resolution, and production settings

[Documentation home](index.md) · [Initial conditions](initial-conditions.md) · [Outputs](outputs.md)

A successful small run checks the workflow, not scientific accuracy. For a production simulation, establish the scales and statistics you need, then vary particle count, force mesh, time stepping, and screening settings until the results converge within your chosen tolerance.

## Reference production choices

The original manual suggested the following starting points for its tested COLA models:

| Setting | Example starting value | Interpretation |
| --- | --- | --- |
| `simulation_boxsize` | `512.0` or `1024.0` | Box length in `Mpc/h`, not volume. |
| `particle_Npart_1D` | `512` or `1024`, respectively | Number along each axis; total is its cube. |
| `force_nmesh` | `3 * particle_Npart_1D` | Three force cells per initial interparticle separation. |
| `timestep_nsteps` | `{50}` | Example total number of integration steps. |
| `timestep_scalefactor_spacing` | `"linear"` | Steps spaced linearly in scale factor. |

These correspond to about one particle per `(Mpc/h)^3` initially. They are expensive production examples, not laptop defaults. In particular, `particle_Npart_1D` is **not** the total particle count; do not cube a desired per-axis count before assigning it. The [backend guide](backend.md) provides small-grid settings for checking execution.

The mesh-to-particle ratio and linear scale-factor spacing follow the COLA accuracy study cited in [references](references.md). They are starting choices, not universal optimal parameters for every gravity model or statistic.

## Spatial and mass resolution

Let $L$ be the box length and $N_{\mathrm{1D}}$ the particle count per axis. Then

$$
N_{\mathrm{tot}}=N_{\mathrm{1D}}^3,\qquad
\ell_{\mathrm{part}}=\frac{L}{N_{\mathrm{1D}}},\qquad
\Delta x_{\mathrm{force}}=\frac{L}{N_{\mathrm{force}}}.
$$

For equal-mass matter particles with the simulated density $\Omega_{\mathrm{particles},0}$,

$$
m_{\mathrm{part}}=\Omega_{\mathrm{particles},0}\rho_{\mathrm{crit},0}
\left(\frac{L}{N_{\mathrm{1D}}}\right)^3.
$$

For the supported frontend setup without massive neutrinos, $\Omega_{\mathrm{particles},0}=\Omega_{m0}$. Using $\rho_{\mathrm{crit},0}\simeq2.775\times10^{11}h^2\,M_\odot\,\mathrm{Mpc}^{-3}$ and $L$ in `Mpc/h`,

$$
m_{\mathrm{part}}\simeq8.325\times10^{10}
\left(\frac{\Omega_{m0}}{0.3}\right)
\left(\frac{L}{1000\,h^{-1}\mathrm{Mpc}}\right)^3
\left(\frac{1000}{N_{\mathrm{1D}}}\right)^3 h^{-1}M_\odot.
$$

Increasing particle count improves mass sampling; increasing the mesh improves force resolution. They address different errors. Mesh and particle memory both grow cubically with their respective per-axis counts, and the solver maintains multiple arrays/grids.

## What limits trustworthy scales?

Three effects highlighted by the manual are:

1. **Time resolution:** too few steps degrade small-scale evolution. More steps improve results until other errors dominate.
2. **Force resolution:** the force-cell size $L/N_{\mathrm{force}}$ limits the resolved density and gravitational force. Changing the box at fixed mesh changes this scale.
3. **Modified-gravity approximations:** the approximate screened force and LPT treatment have model-dependent errors. Increasing mesh resolution alone cannot remove those approximations.

Box size additionally affects the available large-scale modes and finite-volume variance. Particle sampling and shot noise matter at small scales. Determine the usable $k$ range by comparisons and convergence tests at the required accuracy, rather than choosing a single universally trusted cutoff.

For matched GR/modified-gravity ratios, use the same initial phases, seed, normalisation choices, box, particle count, meshes, and output redshifts. For absolute spectra, also validate against an independent calculation or appropriate high-resolution reference. Preserve the parameters, input files, log, source revision, and environment for each comparison.

## When to save particle snapshots

If only the backend's on-the-fly matter spectra are needed, `output_particles = false` can substantially reduce storage. Save particles when you need to reuse exact initial conditions, visualise the density field, or perform post-processing such as an independent halo catalogue. Configure `output_redshifts`, `output_folder`, and `output_fileformat` and verify that the files exist. On-the-fly halo finding can be used separately; saving snapshots is necessary for later external analyses, not for every backend statistic.
