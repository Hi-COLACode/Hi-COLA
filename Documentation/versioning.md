# Versioning and stability

Hi-COLA follows [Semantic Versioning](https://semver.org/) (`MAJOR.MINOR.PATCH`):

- **MAJOR**: breaking changes to the public Python API or to parameter-file
  formats that existing runs rely on.
- **MINOR**: new features added in a backwards-compatible way.
- **PATCH**: bug fixes and documentation corrections with no behavioural change.

Until `1.0.0` is tagged, Hi-COLA is pre-release: interfaces may still change
between minor versions, and the guarantees below apply to the supported surface
rather than every internal module.

## Stability tiers

| Tier | Contents | Guarantee |
|------|----------|-----------|
| **Stable** | Python frontend (`HiCOLA.Frontend`), the `hicola-frontend` pip package, the COLASolver backend in the `hi-cola` conda packages and native builds (`HICOLA_BUILD_NATIVE=1`), `HiCOLA.backend`, and the frontend and backend parameter-file formats | Breaking changes only in a MAJOR release. |
| **Reference / legacy** | `TestData` submodule, bundled example datasets | Not a runtime dependency. Provided for reference and regression work only. |

## Submodules

Installing Hi-COLA does not need `--recurse-submodules`. The `TestData`
submodule is reference data and is not used by the Python frontend, the native
build or the standard test suite.
