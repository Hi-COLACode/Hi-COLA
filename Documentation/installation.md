# Installation

[Documentation home](index.md) · [Next: quickstart](quickstart.md)

## Install a released package

Hi-COLA is distributed under two package names. Both are imported in Python as `HiCOLA`.

| Method | Package | Contents |
| --- | --- | --- |
| pip | `hicola-frontend` | Python frontend and plotting tools. No N-body solver. |
| conda | `hi-cola` | Python frontend, plotting tools, and the compiled COLASolver executable. |

### pip

```sh
python -m pip install hicola-frontend
```

The wheel installs the frontend, its Python dependencies (NumPy, SciPy, SymPy, Matplotlib, ConfigObj, and pandas), and the [`hicola_plot_*` commands](plotting.md). Python 3.9 to 3.12 is supported. The wheel does not contain the compiled backend; use conda or a source build to run simulations.

### conda

```sh
conda install -c conda-forge hi-cola
```

The conda package additionally contains the serial COLASolver executable, so frontend tables can be generated and simulations run [from Python](backend.md#run-from-python) without compiling C++ code. Packages are built for Linux, macOS, and Windows.

Development builds from the default branch are uploaded to a test label rather than conda-forge:

```sh
conda install -c ACCOUNT/label/test -c conda-forge hi-cola
```

Replace `ACCOUNT` with the Anaconda.org account configured for the repository. [RELEASING.md](repo:RELEASING.md) describes version numbering and how releases reach PyPI, TestPyPI, and conda-forge.

### Check the installation

```sh
python -c "import HiCOLA; from HiCOLA import backend; print(backend.is_available())"
hicola_plot_pofk --help
```

`backend.is_available()` prints `True` only when a solver executable is installed (conda or native source build) or selected with `HICOLA_BACKEND_EXECUTABLE`.

## Install from source

Clone the repository to modify the code, use unreleased changes, or build the backend with MPI:

```sh
git clone https://github.com/Hi-COLACode/Hi-COLA.git
cd Hi-COLA
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

An SSH checkout can use `git@github.com:Hi-COLACode/Hi-COLA.git` instead. `pyproject.toml` declares the Python dependencies, so `pip` installs them. Add `'.[tests]'` instead of `.` to also install the test tools. Record your environment for reproducibility:

```sh
python -m pip freeze > environment.txt
```

### Native build through pip

Setting `HICOLA_BUILD_NATIVE=1` makes a source install compile COLASolver with CMake and install it inside the package, as the conda package does:

```sh
HICOLA_BUILD_NATIVE=1 python -m pip install .
```

This requires CMake, a C++17 compiler, FFTW, GSL, and Lua development files. Native builds are serial by default. The CMake option `HICOLA_USE_MPI` additionally requires MPI and MPI-enabled FFTW. A native build cannot be installed in editable mode.

### Makefile build

The backend can also be built directly from `HiCOLA/FML/COLASolver`. The Makefile enables MPI, OpenMP, and FFTW threads by default. With MPI enabled, provide the MPI compiler wrapper and MPI-enabled FFTW library. On a cluster, load the relevant compiler, FFTW, GSL, Lua, and MPI modules first.

Edit the [Makefile](repo:HiCOLA/FML/COLASolver/Makefile) to select compatible compilers, feature flags, and dependency paths. `FML_INCLUDE` must point to the directory **containing** `FML`, which is `HiCOLA` for the bundled copy. The default points to `~/local/FML` and usually needs overriding. Build from the repository root:

```sh
make -C HiCOLA/FML/COLASolver FML_INCLUDE="$PWD/HiCOLA" -j2
```

For a serial build, also pass `USE_MPI=false`. If OpenMP is unavailable, pass `USE_OMP=false`; keep required FFTW, GSL, and Lua support enabled. After changing compiler or feature flags, run `make -C HiCOLA/FML/COLASolver clean` before rebuilding to avoid reusing incompatible object files.

Success produces `HiCOLA/FML/COLASolver/nbody`. To launch it through the Python API, set `HICOLA_BACKEND_EXECUTABLE` to its path. Compilation alone does not validate a simulation: follow the [backend guide](backend.md) with small grids and inspect the output.

### Docker runtime

The repository `Dockerfile` provides an environment with the compilers and libraries needed by the [small example simulations](examples.md).

## Tests

See [testing](testing.md) for the unit-test suite and regression checks, and the [publishing guide](publishing.md) for documentation builds.

## Historical dependency reference

The original quickstart recorded successful use of NumPy 1.23.5, SciPy 1.10.0, SymPy 1.11.1, ConfigObj 5.0.8, FFTW 3.3.8, GSL 2.4, and Lua 5.4.1. This is a historical compatibility record, not a requirement to downgrade to those releases. The documentation build environment is pinned separately in `requirements-docs.lock`; a documentation build is not a scientific backend-validation run.

For Apple Clang, the old guide mentioned a repeated template-default-argument error in `FriendsOfFriends/FoF.h`. If that diagnostic occurs, inspect the declarations for repeated defaults or use a compatible GCC build. Do not apply historical line-number edits blindly to a different revision. FFTW thread library choices must match your installed library and compiler configuration; retain consistent MPI/OpenMP flags across all dependencies.
