"""Launch COLASolver using a Lua file referencing the frontend output tables."""

import logging
import os
from pathlib import Path
import shutil
import subprocess

logger = logging.getLogger(__name__)


EXECUTABLE_NAME = "hi-cola-nbody"


def executable_path():
    """Locate HICOLA_BACKEND_EXECUTABLE, or hi-cola-nbody on PATH."""
    override = os.environ.get("HICOLA_BACKEND_EXECUTABLE")
    if override:
        path = Path(override).expanduser().resolve()
    else:
        found = shutil.which(EXECUTABLE_NAME)
        path = Path(found).resolve() if found else None
    logger.debug("COLASolver executable: %s", path)
    if path is None or not path.is_file():
        raise FileNotFoundError(
            "COLASolver is not installed. Install the hi-cola conda package, put "
            f"{EXECUTABLE_NAME} on PATH, or set HICOLA_BACKEND_EXECUTABLE to a "
            "compiled solver. The Python frontend remains available without COLASolver."
        )
    return path


def is_available():
    """Return whether a solver executable is present."""
    try:
        executable_path()
    except FileNotFoundError:
        return False
    return True


def run(parameter_file, *, cwd=None, capture_output=False, timeout=None):
    """Run the simulation; paths in the Lua file are relative to cwd.

    The parameter file must reference the frontend's expansion and force
    tables. A failed solver process raises subprocess.CalledProcessError.
    """
    parameters = Path(parameter_file).expanduser().resolve()
    if not parameters.is_file():
        raise FileNotFoundError(parameters)
    command = [str(executable_path()), str(parameters)]
    logger.debug("Running COLASolver: %s (cwd=%s)", command, cwd)
    return subprocess.run(
        command, cwd=cwd, check=True, text=True,
        capture_output=capture_output, timeout=timeout,
    )
