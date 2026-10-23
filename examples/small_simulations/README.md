# Small simulations

`run_simulations.py` builds Hi-COLA, prepares model tables and runs GR,
Horndeski, K-mouflage and lightcone examples. Each simulation uses 16³ particles,
a 256 Mpc/h box and five timesteps from redshift 49 to zero. Outputs include
power spectra and particle snapshots; the lightcone example also writes ASCII
particle crossings. These small runs demonstrate execution, not convergence.

The K-mouflage example uses an analytic LCDM expansion and manufactured coupling
and screening tables. It demonstrates the backend, not a calibrated physical
K-mouflage model. The Horndeski example runs the cubic-Galileon frontend.

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

`--work` must be a new directory. Use `--cases gr` to run one example or
`--ranks 2` for two MPI processes. Build commands, process IDs, logs and simulation
outputs are written under `--work`. Compilation writes build products into the
source checkout.

An existing installation can run the same Python command with its dependency
prefix and Lua flags. It needs a C++ compiler, make, MPI, MPI-enabled FFTW, GSL,
Lua and the Python packages in the root `requirements-test.txt`.

`backend_minimal.lua` contains shared settings. The runner appends particle/grid
sizes, file paths and model options. `horndeski.ini` and `numerical.ini` configure
the frontend tables. No comparison code is required to run these examples.
