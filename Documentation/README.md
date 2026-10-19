# Hi-COLA documentation

Hi-COLA connects a Python solver for reduced Horndeski models to a C++ COLA N-body simulator. The frontend exports background and fifth-force tables; the backend uses those tables with an initial matter spectrum to evolve particles.

## Start here

1. [Install the frontend and backend](installation.md).
2. [Generate tables with a supplied model](quickstart.md).
3. [Configure your theory and cosmology](frontend.md).
4. [Run a small simulation, then scale up](backend.md).
5. [Inspect outputs](outputs.md) and [troubleshoot problems](troubleshooting.md).

These pages describe the checked-in implementation. Paths in shell examples are relative to the repository root unless stated otherwise. The quickstart and user-manual material is integrated into these pages, including the theory, equation conventions, initial conditions, accuracy guidance, and output interpretation. Where template comments disagree with executable code, these pages follow the code.

## Reference and development

- [FML COLA solver documentation](repo:HiCOLA/FML/COLASolver/README.md)
- [Frontend model examples](repo:HiCOLA/Parameter_files/Frontend_parameter_files/examples)
- [Backend parameter template](repo:HiCOLA/Parameter_files/Backend_parameter_files/Hi-COLA_params.lua)

For papers, citation guidance, licence, and team contact details, see [citation and references](references.md) and the [project README](repo:README.md). Keep model inputs, generated tables, backend parameters, logs, and the Git revision together so a simulation can be reproduced.

```{toctree}
:maxdepth: 2
:caption: Documentation

installation
quickstart
frontend
theory
equations
backend
initial-conditions
resolution
outputs
troubleshooting
references
api/index
publishing
```
