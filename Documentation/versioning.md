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
| **Reference / legacy** | Bundled example datasets and the FML library's own example programs | Not a runtime dependency. Provided for reference and regression work only. |

## FML example data

Hi-COLA has no git submodules. The example programs under `HiCOLA/FML/*/example`
are part of the bundled FML library and read sample data from `HiCOLA/TestData`,
which is not included. To run them, clone
[HAWinther/TestData](https://github.com/HAWinther/TestData) into `HiCOLA/TestData`.
The frontend, COLASolver and the test suites do not use this data.
