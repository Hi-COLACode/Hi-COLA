"""Launch COLASolver using a Lua file referencing the frontend output tables."""

import logging
import os
from pathlib import Path
import subprocess

logger = logging.getLogger(__name__)


def executable_path():
    """Locate the bundled solver or an explicit HICOLA_BACKEND_EXECUTABLE."""
    override = os.environ.get("HICOLA_BACKEND_EXECUTABLE")
    name = "hi-cola-nbody.exe" if os.name == "nt" else "hi-cola-nbody"
    path = Path(override).expanduser().resolve() if override else Path(__file__).parent / "bin" / name
    if not path.is_file():
        raise FileNotFoundError(
            "COLASolver is not installed. Install a native conda build or set "
            "HICOLA_BACKEND_EXECUTABLE to your compiled nbody executable. "
            "The Python frontend remains available without COLASolver."
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
