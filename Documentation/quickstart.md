# Quickstart: generate simulation inputs

[Documentation home](index.md) · [Installation](installation.md) · [Next: backend](backend.md)

Use a supplied cubic Galileon model before writing a new theory. The blank `horndeski_parameters.ini` is a template and cannot run until `K`, `G3`, and `G4` are filled in.

Activate your [installed environment](installation.md), then run from the repository root:

```sh
MPLBACKEND=Agg python -m HiCOLA.Frontend.generate_simulation_input \
  HiCOLA/Parameter_files/Frontend_parameter_files/examples/cubic_Galileon/cuGal_f0p2.ini \
  HiCOLA/Parameter_files/Frontend_parameter_files/examples/WMAP_cosmological_parameters.ini
```

The command takes exactly two positional arguments: the theory INI followed by the numerical/cosmology INI. This example sets `output_directory = './Hi-COLA_Output'`, resolved against the working directory. It prints the model, cosmology, solved closure parameter, and output directory, then writes:

```text
Hi-COLA_Output/cuGal_f0p2_14dp_WMAP_expansion.txt
Hi-COLA_Output/cuGal_f0p2_14dp_WMAP_force.txt
```

Existing filenames trigger renaming, so use the filenames reported by the run. Both files should have three finite numeric columns, increasing scale factor, and a final scale factor of one. See [table formats and checks](outputs.md).

To customise the example, copy both INI files into a working directory, change `model_name`, `cosmo_name`, and `output_directory`, then edit model and cosmology values using the [frontend reference](frontend.md). Keep the frontend cosmology consistent with backend Lua values.

The generated tables are inputs to the simulator, not particle snapshots. Continue with the [backend guide](backend.md) to configure the Lua file and run a small N-body job.

For theory conventions, see [theory and code structure](theory.md) and [background equations](equations.md). For production choices, follow the [initial-condition](initial-conditions.md) and [resolution](resolution.md) guides after the small backend run succeeds.
