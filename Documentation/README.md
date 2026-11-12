# Hi-COLA documentation

Hi-COLA runs fast, approximate N-body simulations of structure formation in reduced Horndeski gravity. A Python frontend turns a Lagrangian into expansion and fifth-force tables; a C++ COLA backend uses those tables to evolve particles.

1. [Install](installation.md) with `pip` or `conda`.
2. Follow the [quickstart](quickstart.md): generate tables, run small simulations, and plot them.
3. Configure your own [model](frontend.md) and [simulation](backend.md).
4. Interpret the [outputs](outputs.md), or check [troubleshooting](troubleshooting.md).

The [theory page](theory.md) gives the equations solved, and [references](references.md) gives the citation and licence. These pages describe the checked-in code; where template comments disagree with the code, they follow the code.

```{toctree}
:maxdepth: 2
:hidden:

installation
quickstart
frontend
backend
outputs
theory
troubleshooting
references
api/index
development
```
