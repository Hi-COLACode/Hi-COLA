# Hi-COLA

<p align="left">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="Hi-COLA-dark.png">
    <source media="(prefers-color-scheme: light)" srcset="Hi-COLA-light.png">
    <img alt="Hi-COLA logo" src="Hi-COLA-light.png" width="120">
  </picture>
</p>

Hi-COLA is a package that runs fast, approximate N-body simulations of non-linear structure formation in reduced Horndeski gravity (Horndeski theories with luminal gravitational waves). 

Hi-COLA is not hard-coded to solve specific Horndeski theories, but is designed to be generic with respect to the reduced Horndeski class. Given an input Lagrangian, Hi-COLA's frontend dynamically constructs the appropriate field equations and consistently solves for the cosmological background, linear growth, and screened fifth force of that theory. This is passed to the backend, an adaptation of the FML library, where a hybrid N-body simulation at significantly reduced computational and temporal cost compared to traditional N-body codes is run. By analysing the particle snapshots, one is able to study the formation of structure through statistics like the matter power spectrum.

See the Documentation folder for guidance on how to install and use Hi-COLA. There you will find both a quickstart guide and a more detailed user manual.

### Python access to COLASolver

The Python package is imported as `HiCOLA`. Standard pip builds contain the
frontend; native builds also install the COLASolver executable. After generating
the expansion and force tables, reference them in a backend Lua parameter file
and launch the simulation from Python:

```python
from HiCOLA import backend

backend.run("/path/to/run.lua", cwd="/path/to/simulation")
```

Relative paths inside the Lua file are resolved from `cwd`. Solver errors raise
`subprocess.CalledProcessError`. `backend.is_available()` reports whether an
executable is present. A separately compiled solver can be selected with
`HICOLA_BACKEND_EXECUTABLE=/path/to/nbody`.

Native source builds require CMake, a C++17 compiler, FFTW, GSL and Lua. Set
`HICOLA_BUILD_NATIVE=1` when installing from source to compile the solver. Native
builds are serial by default; CMake's `HICOLA_USE_MPI` option additionally requires
MPI and MPI-enabled FFTW. The Python frontend does not implement an N-body solver
to substitute for a missing executable.
EuclidEmulator2 reference predictions are available only when its optional FML
header is present; the bundled sources do not include that emulator.

### Package channels

The build workflow publishes default-branch development versions to TestPyPI.
Release tags use `vMAJOR.MINOR.PATCH`, starting at `v0.0.0`; package managers
normalize the version to `0.0.0`, with `.devN` appended for development builds.
Version-tag pushes publish verified Python artifacts to PyPI using Trusted Publishing.

Native conda builds use the `.github/conda` recipe. In the GitHub repository,
set the `ANACONDA_OWNER` repository Actions variable and create `conda-test` and
`conda` environments for deployment rules. Store `ANACONDA_API_TOKEN`, with upload
access to that account, only as a repository Actions secret under Settings →
Secrets and variables → Actions → Secrets. Do not define the same secret in either
environment: environment secrets override repository secrets. When the token
expires, create a replacement in Anaconda.org and update this one repository
secret; both environments use it.
The workflow builds and tests on Linux, macOS and Windows before uploading each
platform's package to the account's `test` label on default-branch pushes and
`main` label on version-tag pushes. Without `ANACONDA_OWNER`, these jobs are
skipped. Development packages can be installed with:

```sh
conda install -c YOUR_ACCOUNT/label/test -c conda-forge hi-cola
```

Stable releases use `conda install -c YOUR_ACCOUNT -c conda-forge hi-cola`.

This account channel is separate from official conda-forge distribution, which
requires a reviewed recipe and a conda-forge feedstock.

Also see [our paper introducing Hi-COLA](https://iopscience.iop.org/article/10.1088/1475-7516/2023/03/040) [(pre-print version)](https://arxiv.org/abs/2209.01666), where we detail the work gone into creating the first version of Hi-COLA. We encourage the community to use Hi-COLA for their own research, and request that you cite us as detailed in the user manual (a handy bibtex is also provided there). Hi-COLA is provided under a [CC BY 4.0 licence](https://creativecommons.org/licenses/by/4.0/).

Hi-COLA is set up to natively operate with Vainshtein screening. However, with a few specific settings, K-mouflage models can also be simulated -- see our [paper on K-mouflage in Hi-COLA](https://iopscience.iop.org/article/10.1088/1475-7516/2024/11/052) ([pre-print here](https://arxiv.org/abs/2407.00855)). 

How does the accuracy of Hi-COLA fare relative to other approximate tools for LSS? You can find a detailed comparison exercise [in this paper](https://academic.oup.com/mnras/article/536/1/664/7903359?login=true) ([pre-print version here](https://arxiv.org/abs/2406.13667)). 

Hi-COLA remains under active development. You can contact the Hi-COLA team at team.hicola@gmail.com, or at our individual emails (see papers linked above).
