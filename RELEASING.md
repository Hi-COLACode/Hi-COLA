# Release and publishing setup

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
