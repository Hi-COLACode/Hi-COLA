# Installation

[Documentation home](index.md) · [Next: quickstart](quickstart.md)

## Obtain the source

```sh
git clone https://github.com/Hi-COLACode/Hi-COLA.git
cd Hi-COLA
```

An SSH checkout can use `git@github.com:Hi-COLACode/Hi-COLA.git` instead. On a cluster, load the relevant compiler, FFTW, GSL, Lua, and MPI modules before configuring or running the backend.

## Python frontend

From a checkout of Hi-COLA, create an isolated Python 3 environment:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip install -e .
python -c "import numpy, scipy, sympy, matplotlib, configobj, HiCOLA"
```

The frontend imports NumPy, SciPy, SymPy, Matplotlib, and ConfigObj. `setup.py` registers the package but does not declare these dependencies, so install them explicitly. `requirements.txt` lists the frontend runtime imports; it does not pin versions. Record your environment for reproducibility:

```sh
python -m pip freeze > environment.txt
```

For a headless session, prefix frontend commands with `MPLBACKEND=Agg`.

## C++ backend

The backend lives in `HiCOLA/FML/COLASolver`. It requires a C++17 compiler, FFTW, GSL, and Lua development headers and libraries. The default Makefile additionally enables MPI, OpenMP, and FFTW threads. With MPI enabled, provide the MPI compiler wrapper and MPI-enabled FFTW library. Install these dependencies through your machine's package manager or cluster modules; their include and library paths depend on your environment.

Edit the [Makefile](repo:HiCOLA/FML/COLASolver/Makefile) to select compatible compiler, feature flags, and dependency paths. `FML_INCLUDE` must point to the directory **containing** `FML`, which is `HiCOLA` for the bundled copy. The default points to `~/local/FML` and usually needs overriding.

After configuring dependency paths, build from the repository root:

```sh
make -C HiCOLA/FML/COLASolver FML_INCLUDE="$PWD/HiCOLA" -j2
```

For a serial build, also pass `USE_MPI=false`. If OpenMP is unavailable, pass `USE_OMP=false`; keep required FFTW, GSL, and Lua support enabled. Use the same feature flags for subsequent builds. After changing compiler or feature flags, run `make -C HiCOLA/FML/COLASolver clean` before rebuilding to avoid reusing incompatible object files.

Success produces `HiCOLA/FML/COLASolver/nbody`. Compilation alone does not validate a simulation: follow the [backend guide](backend.md) with small grids and inspect the output.

## Documentation checks

See the [publishing guide](publishing.md) for the Sphinx build and documentation-specific tests. These checks are independent of the scientific unit-test branch.

## Historical dependency reference

The original quickstart recorded successful use of NumPy 1.23.5, SciPy 1.10.0, SymPy 1.11.1, ConfigObj 5.0.8, FFTW 3.3.8, GSL 2.4, and Lua 5.4.1. This is a historical compatibility record, not a requirement to downgrade to those releases. The documentation CI environment is pinned separately in `requirements-docs.lock`; a documentation build is not a scientific backend-validation run.

For Apple Clang, the old guide mentioned a repeated template-default-argument error in `FriendsOfFriends/FoF.h`. If that diagnostic occurs, inspect the declarations for repeated defaults or use a compatible GCC build. Do not apply historical line-number edits blindly to a different revision. FFTW thread library choices must match your installed library and compiler configuration; retain consistent MPI/OpenMP flags across all dependencies.
