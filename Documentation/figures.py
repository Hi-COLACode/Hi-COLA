"""Generate the example figures shown in the documentation at build time.

Each figure is exposed to the pages as an include fragment in ``_generated``:
either a MyST figure, or a note explaining why the figure is unavailable.

* The frontend-table figure is produced on every build by running the Horndeski
  frontend on a bundled example model and plotting the result with
  ``hicola_plot_frontend``.
* Simulation figures are copied from a completed run of
  ``examples/small_simulations/run_simulations.py`` whose ``--work`` directory
  is given by the ``DOCS_SIMULATION_WORK`` environment variable. The
  Documentation workflow performs that run before building.
"""
import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

logger = logging.getLogger(__name__)

FRONTEND_INPUTS = (
    "HiCOLA/Parameter_files/Frontend_parameter_files/examples/cubic_Galileon/cuGal_f0p2.ini",
    "HiCOLA/Parameter_files/Frontend_parameter_files/examples/WMAP_cosmological_parameters.ini",
)
SIMULATION_FIGURES = ("pofk", "background", "growth")
SIMULATION_NOTE = (
    "This figure is produced when the documentation build runs the "
    "[small example simulations](examples.md); the Documentation workflow does so. "
    "Set `DOCS_SIMULATION_WORK` to the `--work` directory of a completed run "
    "to include it in a local build."
)


def fragment(directory, name, image, alt, caption, note):
    """Write the include fragment for one figure."""
    if image is not None:
        text = f"```{{figure}} /_generated/{image.name}\n:alt: {alt}\n:width: 100%\n\n{caption}\n```\n"
    else:
        text = f":::{{note}}\n{note}\n:::\n"
    (directory / f"{name}.md").write_text(text)
    logger.debug("Figure fragment %s uses image %s", name, image)


def frontend_figure(root, directory):
    """Run the frontend on the example model and plot its tables."""
    plotting = root / "HiCOLA/Plotting/frontend.py"
    if not plotting.is_file():
        logger.debug("No frontend plotting tool in this checkout")
        return None
    environment = dict(os.environ, MPLBACKEND="Agg",
                       PYTHONPATH=os.pathsep.join(filter(None, [str(root), os.environ.get("PYTHONPATH")])))
    with tempfile.TemporaryDirectory() as work:
        logger.debug("Generating example frontend tables in %s", work)
        subprocess.run([sys.executable, "-m", "HiCOLA.Frontend.generate_simulation_input",
                        *(str(root / path) for path in FRONTEND_INPUTS)],
                       cwd=work, env=environment, check=True, stdout=subprocess.DEVNULL)
        tables = Path(work) / "Hi-COLA_Output"
        subprocess.run([sys.executable, "-m", "HiCOLA.Plotting.frontend", str(tables),
                        "--label", "cubic Galileon example", "--out", work],
                       cwd=work, env=environment, check=True, stdout=subprocess.DEVNULL)
        image = directory / "frontend.png"
        shutil.copyfile(Path(work) / "frontend.png", image)
    return image


def simulation_figure(directory, name):
    """Copy one cubic Galileon simulation plot from a completed example run."""
    work = os.environ.get("DOCS_SIMULATION_WORK")
    if not work:
        return None
    source = Path(work) / "cubic_galileon" / "plots" / f"{name}.png"
    if not source.is_file():
        raise FileNotFoundError(f"DOCS_SIMULATION_WORK has no {source}")
    image = directory / f"simulation-{name}.png"
    shutil.copyfile(source, image)
    return image


def generate(root, documentation):
    """Create all figure fragments under ``documentation/_generated``."""
    directory = documentation / "_generated"
    shutil.rmtree(directory, ignore_errors=True)
    directory.mkdir()
    fragment(directory, "frontend", frontend_figure(root, directory),
             "Frontend tables plotted against scale factor",
             "`hicola_plot_frontend` applied to the tables generated from the "
             "cubic Galileon example model during this documentation build.",
             "The frontend plotting tool is not present in this checkout.")
    captions = {
        "pofk": ("Matter power spectra with ratio panels",
                 "`hicola_plot_pofk` for the cubic Galileon example with `--reference` set to the GR example: "
                 "the spectra at each output redshift and their ratio with ±1% and ±2.5% bands."),
        "background": ("Expansion history and density fractions",
                       "`hicola_plot_background` for the cubic Galileon example against the GR example."),
        "growth": ("Linear growth and effective gravitational coupling",
                   "`hicola_plot_growth` for the cubic Galileon example against the GR example."),
    }
    for name in SIMULATION_FIGURES:
        alt, caption = captions[name]
        caption += (" Generated during this documentation build from 16³-particle example runs; "
                    "it demonstrates the tools, not converged results.")
        fragment(directory, name, simulation_figure(directory, name), alt, caption, SIMULATION_NOTE)
