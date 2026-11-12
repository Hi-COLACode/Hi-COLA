# Quickstart

## Generate frontend tables

From a source checkout, run the supplied cubic Galileon model:

```sh
python -m HiCOLA.Frontend.generate_simulation_input \
  HiCOLA/Parameter_files/Frontend_parameter_files/examples/cubic_Galileon/cuGal_f0p2.ini \
  HiCOLA/Parameter_files/Frontend_parameter_files/examples/WMAP_cosmological_parameters.ini
```

The arguments are the theory file and the numerical/cosmology file. The run prints the solved closure parameter and writes `cuGal_f0p2_14dp_WMAP_{expansion,force,full}.txt` to `Hi-COLA_Output/` in the working directory; existing names are renamed, so use the printed paths. The [outputs page](outputs.md#frontend-tables) describes the columns:

```{include} _generated/frontend.md
```

To use your own model, copy both INI files, change `model_name`, `cosmo_name`, and `output_directory`, and edit them following the [frontend reference](frontend.md).

## Run the example simulations

`examples/small_simulations` contains complete runs with 16³ particles in an 80 Mpc/h box, a 32³ force mesh, and five steps from redshift 49 to 0:

| Case | Model |
| --- | --- |
| `gr` | GR |
| `cubic_galileon` | Cubic Galileon, tables from the Horndeski frontend |
| `kmouflage` | K-mouflage with manufactured tables (demonstrates the backend, not a calibrated model) |
| `lightcone` | GR lightcone from redshift 0.25 with box replication; the slowest case |

These runs demonstrate execution, not convergence. In Docker, from the repository root:

```sh
docker build -t hicola-runtime .
docker run --rm -it --user "$(id -u):$(id -g)" -v "$PWD:/workspace" hicola-runtime
python examples/small_simulations/run_simulations.py \
  --work /workspace/example-output --prefix /usr \
  --make-arg LUA_INCLUDE=/usr/include/lua5.4 --make-arg 'LUA_LINK=-llua5.4 -ldl'
```

The runner builds COLASolver, generates tables, runs each case, and plots it into `<case>/plots/`, comparing modified-gravity cases with GR. `--work` must be a new directory. `--cases gr cubic_galileon` selects cases and `--ranks 2` uses two MPI processes. Outside Docker it needs a C++ compiler, make, MPI, MPI-enabled FFTW, GSL, Lua, and `python -m pip install -e .`.

```{include} _generated/pofk.md
```

`backend_minimal.lua` holds the shared settings, and `gr.lua`, `cubic_galileon.lua`, `kmouflage.lua`, and `lightcone.lua` select model and outputs. A built solver can run one directly; create its `output_folder` first:

```sh
mkdir -p Hi-COLA_Output/gr/output
mpiexec -n 1 HiCOLA/FML/COLASolver/nbody examples/small_simulations/gr.lua
```

Direct modified-gravity runs read tables from `Hi-COLA_Output/frontend/`. Continue with the [backend guide](backend.md) to set up your own simulation.
