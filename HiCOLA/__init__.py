"""Hi-COLA frontend and optional compiled simulation backend."""

from . import backend

USING_CPP = backend.is_available()

__all__ = ["backend", "USING_CPP"]
