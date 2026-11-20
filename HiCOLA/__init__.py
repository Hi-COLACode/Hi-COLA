"""Hi-COLA frontend and optional compiled simulation backend."""

from importlib.metadata import PackageNotFoundError, version

from . import backend

try:
    __version__ = version("hicola-frontend")
except PackageNotFoundError:
    # Imported from a source tree that has not been installed.
    __version__ = "0.0.0"

USING_CPP = backend.is_available()

__all__ = ["backend", "USING_CPP", "__version__"]
