"""Check frontend process launching and missing-backend behaviour."""

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from HiCOLA import backend


class BackendTests(unittest.TestCase):
    def test_missing_backend(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch.dict(os.environ, {"HICOLA_BACKEND_EXECUTABLE": str(Path(directory) / "missing")}):
                self.assertFalse(backend.is_available())
                with self.assertRaisesRegex(FileNotFoundError, "COLASolver is not installed"):
                    backend.executable_path()

    def test_launch_uses_parameter_file_and_working_directory(self):
        with tempfile.TemporaryDirectory(prefix="hicola space ") as directory:
            script = Path(directory) / "parameter file.py"
            script.write_text("from pathlib import Path\nprint(Path.cwd())\n", encoding="utf-8")
            with patch.dict(os.environ, {"HICOLA_BACKEND_EXECUTABLE": sys.executable}):
                result = backend.run(script, cwd=directory, capture_output=True)
            self.assertEqual(Path(result.stdout.strip()).resolve(), Path(directory).resolve())

    def test_failed_process_is_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            script = Path(directory) / "failure.py"
            script.write_text("raise SystemExit(7)\n", encoding="utf-8")
            with patch.dict(os.environ, {"HICOLA_BACKEND_EXECUTABLE": sys.executable}):
                with self.assertRaises(subprocess.CalledProcessError) as caught:
                    backend.run(script, capture_output=True)
            self.assertEqual(caught.exception.returncode, 7)

    def test_missing_parameter_file(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(FileNotFoundError):
                backend.run(Path(directory) / "missing.lua")
