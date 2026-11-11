# Small simulations

`run_simulations.py` builds Hi-COLA, prepares model tables and runs GR,
cubic Galileon, K-mouflage and lightcone examples. Each simulation uses 16³ particles,
a box of 80 Mpc/h, a 32³ force mesh and five timesteps from redshift 49 to zero. Outputs include
power spectra and particle snapshots; the lightcone example also writes ASCII
particle crossings. These small runs demonstrate execution, not convergence.
The lightcone starts at redshift 0.25 with box replication enabled to cover its
radial extent beyond the original box. Replication makes this example more
expensive than the other small simulations.

The K-mouflage example uses an analytic LCDM expansion and manufactured coupling
and screening tables. It demonstrates the backend, not a calibrated physical
K-mouflage model. The cubic Galileon example generates its tables with the Horndeski frontend.

Build and enter the runtime environment from the repository root:

```sh
docker build -t hicola-runtime .
docker run --rm -it --user "$(id -u):$(id -g)" \
  -v "$PWD:/workspace" hicola-runtime
```

Inside the container, run all four examples:

```sh
python examples/small_simulations/run_simulations.py \
  --work /workspace/example-output --prefix /usr \
  --make-arg LUA_INCLUDE=/usr/include/lua5.4 \
  --make-arg 'LUA_LINK=-llua5.4 -ldl'
```

After each simulation the script draws its plots into `<case>/plots/` with the
`hicola_plot_*` tools installed with the package:

| Tool | Reads | Draws |
|---|---|---|
| `hicola_plot_pofk` | `pofk*_z*.txt` | P(k) at every output redshift |
| `hicola_plot_background` | `cosmology_*.txt` | H/H0, d ln H / d ln a and density fractions against a |
| `hicola_plot_growth` | `gravitymodel_*_k*.txt` | D1(a, k) and G_eff/G against a for each k |
| `hicola_plot_frontend` | frontend `*_expansion.txt`, `*_force.txt` | E, E'/E, screening coefficient, coupling and G_G4/G_N against a |

Modified-gravity cases are plotted with `--reference` set to the GR run. With a
reference, each panel overlays the run (solid blue) and the reference (dashed orange),
P(k) and growth get one panel per redshift or k, and ratio panels (the boost P/P_GR,
H/H_GR, D1/D1_GR) show ±1% and ±2.5% bands. `--redshifts` and `--wavenumbers` select
which redshifts or k values to plot. Each tool also runs on its own, for example:

```sh
hicola_plot_pofk example-output/cubic_galileon/output \
  --reference example-output/gr/output --bands 0.01 0.025
```

`--bands` takes at most two fractions; give it no values to draw none. Run any tool
with `--help` for its options. Frontend plots are written to
`frontend/<model>/plots/`.

`--work` must be a new directory. Use `--cases gr` to run one example or
`--ranks 2` for two MPI processes. Build commands, process IDs, logs and simulation
outputs are written under `--work`. Compilation writes build products into the
source checkout.

An existing installation can run the same Python command with its dependency
prefix and Lua flags. It needs a C++ compiler, make, MPI, MPI-enabled FFTW, GSL,
Lua and the Python dependencies declared in `pyproject.toml` (`python -m pip install -e .`).

`backend_minimal.lua` contains the shared simulation settings, including box size,
particle count, meshes and timesteps. `gr.lua`, `cubic_galileon.lua`,
`kmouflage.lua` and `lightcone.lua` load it with `dofile` and select their model
and output settings. The Python runner builds the code, prepares frontend tables
and launches the simulations; its generated Lua wrapper supplies only file paths.
`cubic_galileon.ini` and `numerical.ini` configure the frontend tables.

A built backend can also load an example directly from the repository root:

```sh
mkdir -p Hi-COLA_Output/gr/output
mpiexec -n 1 HiCOLA/FML/COLASolver/nbody examples/small_simulations/gr.lua
```

The modified-gravity examples require the expansion and force tables named in
their Lua files. The runner prepares these under `<work>/frontend/` and supplies
their paths automatically; direct runs use `Hi-COLA_Output/frontend/` by default.
Create the output directory specified by `output_folder` before a direct run.
No comparison code is required to run these examples.

The cubic Galileon frontend accepts the example INI files directly:

```sh
python -m HiCOLA.Frontend.generate_simulation_input \
  examples/small_simulations/cubic_galileon.ini \
  examples/small_simulations/numerical.ini
```

`output_directory` in the model INI selects where the tables are written. The
runner supplies its own output directory. Numerical settings come from
`numerical.ini`. K-mouflage uses the runner's manufactured-table generator rather
than the Horndeski frontend.
