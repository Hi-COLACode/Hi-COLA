# Release and publishing setup

The pip distribution is `hicola-frontend` (Python import: `HiCOLA`). Its wheel
contains the frontend; its source archive also contains the C++ backend sources.
The native conda package is named `hi-cola`.

The version comes from git tags through setuptools-scm; no file holds it.
To release, push a `vMAJOR.MINOR.PATCH` tag (for example `git tag v2.0.0 &&
git push upstream v2.0.0`). Packaged releases start at `v2.0.0`, following the
2023 release `Hi-COLA_v1.0.0` (tag `first_release`). Version-tag pushes publish
verified Python artifacts to PyPI using Trusted Publishing.

Untagged commits build as development versions: the next patch after the latest
release tag with `.devN`, where N counts commits since that tag (for example
`2.0.1.dev3`). Default-branch pushes publish them to TestPyPI. Nothing needs
changing after a release.

Before the first release, development versions count from a `v2.0.0.dev0` tag
on the default branch, giving `2.0.0.devN`. Push it once:
`git tag v2.0.0.dev0 <commit> && git push upstream v2.0.0.dev0`. Tags of the
form `vX.Y.Z.devN` set development numbering only and do not publish anything.

Register Trusted Publishers separately on PyPI and TestPyPI for project
`hicola-frontend`, GitHub owner `Hi-COLACode`, repository `Hi-COLA`, and workflow
`build.yml`. Use the `pypi` environment on PyPI and `testpypi` on TestPyPI.
A publisher registered for `Hi_COLA` does not cover this project name.

Native conda builds use the `.github/conda` recipe. Its run requirements repeat the
Python dependencies in `pyproject.toml` with the same minimum versions; the build
workflow checks they match with `.github/conda/check_dependencies.py`. In the GitHub repository,
set the `ANACONDA_OWNER` repository Actions variable and create the `conda-test`
environment for deployment rules. Store `ANACONDA_API_TOKEN`, with upload
access to that account, only as a repository Actions secret under Settings →
Secrets and variables → Actions → Secrets. Do not define the same secret in that
environment: environment secrets override repository secrets. When the token
expires, create a replacement in Anaconda.org and update this one repository
secret; the development publishing job uses it.
The workflow builds and tests the conda package on Linux and macOS for pull
requests, manual runs and default-branch pushes. Default-branch pushes then
upload each platform's package to the account's `test` label. Version tags do
not build or upload conda packages in this workflow. Without `ANACONDA_OWNER`,
the upload is skipped. Development packages can be installed with:

```sh
conda install -c YOUR_ACCOUNT/label/test -c conda-forge hi-cola
```

## Submit a stable release to conda-forge

After this workflow is merged and a tagged release has published to PyPI:

1. Confirm the release's source archive (`.tar.gz`) is available on the
   `hicola-frontend` PyPI project and
   contains `CMakeLists.txt`, the FML C++ sources and parameter files. The pure
   Python wheel is not the source for the native conda package.
2. Fork [conda-forge/staged-recipes](https://github.com/conda-forge/staged-recipes)
   and create a branch from its current `main`. Create `recipes/hi-cola/` and copy
   `.github/conda/meta.yaml`, `.github/conda/build.sh` and
   `.github/conda/run_test.py` there.
3. Adapt the copied recipe for a standalone build:
   - Replace the `HICOLA_VERSION` environment lookup with the published version
     (for example, `{% set version = "2.0.0" %}`).
   - Replace the source environment lookups with the public PyPI source archive
     URL and its SHA256, available on the release's PyPI download page. Keep the
     version variable in the URL so the update bot can substitute new versions.
   - Copy `.github/conda/build.sh` as well, and keep the C++ compiler, native
     dependencies and installed simulation test.
   - Do not copy `.github/conda/conda_build_config.yaml`; conda-forge's global
     pinning supplies the `stdlib('c')` settings. Do not use `noarch`.
   - Complete `about` with project links, and keep its `license` and
     `license_file` entries. Add the license file of any third-party code
     compiled into the solver.
   - Add `extra: recipe-maintainers:` with the maintainers' GitHub usernames.
4. Follow the current [submission checklist](https://conda-forge.org/docs/maintainer/adding_pkgs/),
   test the recipe locally where practical, then open a pull request yourself
   against `conda-forge/staged-recipes`. Address its CI checks and review comments.
5. Once accepted, conda-forge creates `hi-cola-feedstock` and builds and publishes
   the native packages. Verify installation and the small simulation using:

   ```sh
   conda create -n hicola-release -c conda-forge python=3.11 hi-cola
   conda activate hicola-release
   python -c "import HiCOLA; assert HiCOLA.USING_CPP"
   ```

   Run `.github/conda/run_test.py` from outside the source checkout using that
   environment to verify the installed frontend-to-backend simulation.

Users with conda-forge configured can then run `conda install hi-cola`; others
need `conda install -c conda-forge hi-cola`.

For later releases, push a version tag to publish the new source archive to
PyPI. Conda-forge's bot periodically detects new versions and proposes a recipe
version/checksum update in the feedstock. Review the requirements and native
build results, then merge the update there. Feedstock CI handles publication;
this repository does not need a conda-forge upload token. See the
[feedstock maintenance guide](https://conda-forge.org/docs/maintainer/updating_pkgs/).
