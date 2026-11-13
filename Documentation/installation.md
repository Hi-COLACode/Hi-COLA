# Installation

Hi-COLA is imported in Python as `HiCOLA` and is published under two package names:

| Command | Contents |
| --- | --- |
| `python -m pip install hicola-frontend` | Python frontend and [plotting tools](outputs.md#plotting-tools). No N-body solver. |
| `conda install -c conda-forge hi-cola` | The above plus the compiled COLASolver, ready to [run from Python](backend.md#run-the-simulation). |

pip supports Python 3.9 to 3.12 and installs NumPy, SciPy, SymPy, Matplotlib, ConfigObj, and pandas. Conda packages are built for Linux, macOS, and Windows. Development builds of the default branch are on a test label: `conda install -c ACCOUNT/label/test -c conda-forge hi-cola`, with `ACCOUNT` the repository's Anaconda.org account.

Check the installation:

```sh
python -c "from HiCOLA import backend; print(backend.is_available())"
hicola_plot_pofk --help
```

`is_available()` is `True` only when a solver is installed or selected with `HICOLA_BACKEND_EXECUTABLE`.

## From source

Use a checkout to modify the code or to build the solver with MPI:

```sh
git clone https://github.com/Hi-COLACode/Hi-COLA.git
cd Hi-COLA
python -m pip install -e .            # add [tests] for the test tools
```

**Solver through pip.** `HICOLA_BUILD_NATIVE=1 python -m pip install .` compiles COLASolver with CMake into the package, as conda does. It needs CMake, a C++17 compiler, FFTW, GSL, and Lua development files. Builds are serial unless the CMake option `HICOLA_USE_MPI` is set, which also needs MPI and MPI-enabled FFTW. Native builds cannot be editable.

**Solver through make.** The [Makefile](repo:HiCOLA/FML/COLASolver/Makefile) enables MPI, OpenMP, and FFTW threads by default:

```sh
make -C HiCOLA/FML/COLASolver FML_INCLUDE="$PWD/HiCOLA" -j2
```

`FML_INCLUDE` is the directory containing `FML`. Pass `USE_MPI=false` or `USE_OMP=false` to disable those features, and run `make … clean` after changing flags. Edit the Makefile's dependency paths for your system or cluster modules. The result is `HiCOLA/FML/COLASolver/nbody`; point `HICOLA_BACKEND_EXECUTABLE` at it to use it from Python.

**Docker.** The repository `Dockerfile` provides all build and run dependencies; the [quickstart](quickstart.md#run-the-example-simulations) uses it.

Earlier releases were tested with NumPy 1.23.5, SciPy 1.10.0, SymPy 1.11.1, ConfigObj 5.0.8, FFTW 3.3.8, GSL 2.4, and Lua 5.4.1. With Apple Clang, a repeated template default argument in `FriendsOfFriends/FoF.h` may need removing, or use GCC.
