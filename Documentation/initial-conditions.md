# Initial conditions

[Documentation home](index.md) · [Backend](backend.md) · [Accuracy and resolution](resolution.md)

Hi-COLA can generate initial particles from a linear spectrum or read an existing GADGET snapshot. Match the background cosmology, gravity choice, spectrum epoch, and growth normalisation before comparing simulations.

## Spectrum inputs and backscaling

A target-epoch spectrum is scaled to the simulation's initial epoch using the backend growth factors. Supplying a late-time spectrum helps carry information from a linear calculation that includes radiation, which the pressureless particle evolution does not reproduce by itself. Ideally, use a linear Einstein-Boltzmann calculation for the same modified-gravity theory.

For a two-column spectrum, set:

```lua
ic_type_of_input = "powerspectrum"
ic_input_filename = "/absolute/path/model_matterpower.dat"
ic_input_redshift = 0.0
ic_initial_redshift = 49.0
ic_use_gravity_model_GR = false
ic_random_field_type = "gaussian"
```

The current reader expects one header line followed by $k$ in `h/Mpc` and $P(k)$ in `(Mpc/h)^3`. The accepted name is **`powerspectrum`**, not `matterpowerspectrum`. Other supported inputs include `transferfunction` and `transferinfofile`; the latter describes a collection of transfer data rather than a two-column spectrum. Use the [backend input examples](repo:HiCOLA/FML/COLASolver/examples/example_run_different_input) for those formats.

For the power-spectrum input, the implementation uses

$$
P_{\mathrm{ini}}(k)=
\left[\frac{D_1(k,a_{\mathrm{ini}})}{D_1(k,a_{\mathrm{input}})}\right]^2
P_{\mathrm{input}}(k).
$$

`ic_use_gravity_model_GR = false` uses the selected gravity model for this operation. Setting it true uses GR for initial-condition backscaling, which is useful for deliberately matching GR and modified-gravity initial conditions but changes the comparison being made. `ic_sigma8_normalization` can additionally change the input normalisation; record that choice along with the spectrum.

If a theory-specific Einstein-Boltzmann spectrum is unavailable, the original quickstart proposed the approximate construction

$$
P_L^{\mathrm{MG,quasi-back}}(k,z_{\mathrm{ini}})=
\left[\frac{D_1^{\mathrm{MG}}(k,z_{\mathrm{ini}})}
{D_1^{\Lambda\mathrm{CDM}}(k,z_{\mathrm{target}})}\right]^2
P_L^{\Lambda\mathrm{CDM}}(k,z_{\mathrm{target}}).
$$

The two growth functions must use compatible early-time normalisations. This is an approximation, not a replacement for the same theory's linear calculation. If supplying an already backscaled spectrum, set `ic_input_redshift` to its actual epoch so it is not mistakenly backscaled a second time.

## Existing GADGET particles

Snapshot reading is selected by **`ic_random_field_type`** in the current implementation:

```lua
ic_random_field_type = "read_particles"
ic_reconstruct_gadgetfilepath = "/absolute/path/snapshot/gadget"
ic_reconstruct_assigment_method = "CIC"
ic_reconstruct_interlacing = false
ic_reconstruct_smoothing_filter = "sharpk"
```

The GADGET path is the file prefix expected by the reader, conventionally ending in `/gadget`. Use the snapshot's actual initial redshift, units, cosmology, and box size. Keep other required input and growth settings from the Lua template; selecting particle reading does not eliminate every spectrum-related parameter read during initialisation.

When COLA is enabled, the code reconstructs the LPT fields from the supplied particle distribution. Choose `ic_nmesh` consistently with the initial particle grid and inspect the reconstruction settings before comparing to an external simulation. The key `ic_reconstruct_assigment_method` contains that spelling in the actual parameter reader.

## Ordinary Gaussian realisations

To average independent Gaussian runs, vary `ic_random_seed` while keeping the scientific settings fixed:

```lua
ic_random_field_type = "gaussian"
ic_fix_amplitude = false
ic_reverse_phases = false
```

Average the resulting statistics across realisations rather than expecting individual large-scale spectra to match exactly. `ic_random_generator` also affects reproducibility and should be recorded. The backend offers non-Gaussian settings, but they require their own scientific justification and parameter choices.

## Fixed-and-paired runs

Paired-fixed initial conditions reduce variance in estimates of mean large-scale clustering using two simulations. Keep the same seed and all other settings, fix the Fourier amplitudes, then reverse phases in the second run:

| Parameter | First run | Paired run |
| --- | --- | --- |
| `ic_random_seed` | `1234` | `1234` |
| `ic_random_field_type` | `"gaussian"` | `"gaussian"` |
| `ic_fix_amplitude` | `true` | `true` |
| `ic_reverse_phases` | `false` | `true` |

Use different simulation labels and output folders. Average the two spectra to obtain the paired estimate. Pairing does not remove every source of uncertainty, and amplitude-fixed ensembles do not reproduce the ordinary Gaussian ensemble's covariance. Do not use this procedure as a shortcut for estimating covariance matrices. The [references](references.md) include the paired-fixed and covariance studies cited by the original manual.

Save initial snapshots when later reruns must begin from exactly the same particles, especially when compiler, dependency, or random-generator differences could affect regeneration. Set `output_particles = true` and include the desired initial epoch in the output schedule; confirm that the requested snapshot was actually produced.
