# Plotting tools

[Documentation home](index.md) · [Outputs](outputs.md) · [Examples](examples.md)

Installing Hi-COLA with pip or conda adds four command-line tools:

| Command | Reads | Draws |
| --- | --- | --- |
| `hicola_plot_pofk` | Backend `pofk*_z*.txt` | $P(k)$ at every output redshift |
| `hicola_plot_background` | Backend `cosmology_*.txt` | $H/H_0$, $d\ln H/d\ln a$, and density fractions against $a$ |
| `hicola_plot_growth` | Backend `gravitymodel_*_k*.txt` | $D_1(a,k)$ and $G_{\mathrm{eff}}/G$ against $a$ for each $k$ |
| `hicola_plot_frontend` | Frontend `*_expansion.txt` and `*_force.txt` | $E$, $E'/E$, screening coefficient, coupling, and $G_{G_4}/G_N$ against $a$ |

Each command takes the folder holding those files. For example:

```sh
hicola_plot_pofk example-output/cubic_galileon/output \
  --reference example-output/gr/output --bands 0.01 0.025
```

## Options

| Option | Effect |
| --- | --- |
| `--reference FOLDER` | Overlays a reference run of the same kind (dashed) on this run (solid) and adds ratio panels: $P/P_{\mathrm{ref}}$, $H/H_{\mathrm{ref}}$, $D_1/D_{1,\mathrm{ref}}$, or $E/E_{\mathrm{ref}}$. With a reference, $P(k)$ and growth plots get one panel per shared redshift or $k$. |
| `--bands F [F]` | Up to two fractional bands drawn around a ratio of one. The default is `0.01 0.025` (±1% and ±2.5%); give the option no values to draw none. |
| `--redshifts Z ...` | `hicola_plot_pofk` only: select output redshifts. |
| `--wavenumbers K ...` | `hicola_plot_growth` only: select $k$ values. |
| `--label`, `--reference-label` | Names used in titles and legends. Defaults come from the folder names, skipping a trailing `output` folder. |
| `--out FOLDER` | Destination folder. The default is `plots/` beside an `output` folder, otherwise inside the input folder. |
| `--format {png,pdf,svg}` | Figure format; the default is `png`. |

Run any command with `--help` for its complete options. In a headless session, set `MPLBACKEND=Agg`. The [Python API](api/index.rst) lists the `main` functions behind each command.
