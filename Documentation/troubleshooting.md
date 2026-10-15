# Troubleshooting

[Documentation home](index.md) · [Installation](installation.md) · [Frontend](frontend.md) · [Backend](backend.md)

| Symptom | Check and action |
| --- | --- |
| `ModuleNotFoundError` | Activate the intended environment; install `requirements.txt` and the package with that environment's Python. |
| Symbolic parsing fails | Use a complete example, fill `K`, `G3`, `G4`, use explicit multiplication, and declare constants. |
| Wrong parameter substitution | Check that `_symbol` and `_parameter` entries have the same order and count. |
| Closure fails or selects an unexpected branch | Inspect the guess, normalisation, and physical parameter choices; compare to the supplied example and check the printed solved value. |
| Undefined `f_phi` during Hubble closure | Closure index `0` has a known defect; use a supported model with scalar-derivative closure or investigate the solver before proceeding. |
| Index error with closure indices `2` or `3` | The current reader only stores two initial values; template comments overstate supported indices. |
| Direct matter density appears wrong | Keep `Omega_m0 = None` and supply physical baryon/CDM densities; the current direct-input branch reads the radiation key. |
| Solver warnings, NaNs, or infinities | Set `suppression_flag = False`; inspect the first failure, model domain, closure root, and numerical settings. Reject invalid tables. |
| Generated files cannot be found | Resolve `output_directory` against the command's working directory and use the filenames printed after renaming. |
| Renaming raises an exception | Choose a fresh output directory or model name; the helper has limited suffix handling and a 100-iteration limit. |
| `FML/...` header or source not found | Set `FML_INCLUDE` to the directory containing `FML`, normally the checkout's absolute `HiCOLA` path. |
| Missing FFTW/GSL/Lua headers or link symbols | Match development headers, library paths, compiler, and build flags. MPI builds need MPI FFTW as well as MPI itself. |
| Compiler rejects OpenMP | Configure an OpenMP-capable compiler or build with `USE_OMP=false`; clean objects before rebuilding. |
| Force file is missing or its path is empty | Use `HiCOLA_preforce_filename`; replace the obsolete key in the Lua template. |
| Backend reports a missing parameter | With `all_parameters_must_be_in_file = true`, retain keys read for the selected modes; use `ReadParameters.h` as the reference. |
| Run exhausts memory | Begin with the small-grid workflow; particle and mesh counts scale cubically. Account for working grids and particle allocation overhead. |
| Unexpected spectrum or background | Match frontend and backend cosmologies, spectrum epoch, GR back-scaling, sigma8 normalisation, and table coverage. |

For a reproducible report, include the Git revision, operating system/compiler or Python environment, exact command and working directory, both INI files, Lua file, and the first relevant error with nearby log lines. State whether table generation, backend compilation, or simulation execution fails.

Known limitations include Hubble closure argument substitution and the `G3` derivative with respect to `phi`. Explicitly `phi`-dependent models need independent validation of that derivative path.
