"""Build the Python frontend, optionally including the COLASolver executable."""

from email.parser import Parser
import os
from pathlib import Path
import re
import sys

from setuptools import find_packages

native = os.environ.get("HICOLA_BUILD_NATIVE", "0") == "1"
if native:
    from skbuild import setup
else:
    from setuptools import setup


def resolve_version():
    """Return the base version; setup.cfg's tag_build adds any .devN suffix."""
    pkg_info = Path(__file__).resolve().with_name("PKG-INFO")
    if pkg_info.is_file():
        # An sdist records its version, so every build from it agrees.
        metadata = Parser().parsestr(pkg_info.read_text(encoding="utf-8"))
        return metadata["Version"].split(".dev")[0]
    if os.environ.get("GITHUB_REF_TYPE") == "tag":
        tag = os.environ.get("GITHUB_REF_NAME", "")
        if not re.fullmatch(r"v\d+\.\d+\.\d+", tag):
            raise ValueError("Release tags must have the form vMAJOR.MINOR.PATCH")
        return tag[1:]
    return os.environ.get("PKG_VERSION", "0.0.0")


version = resolve_version()

options = {}
if native:
    options["cmake_args"] = [
        "-DHICOLA_BUILD_NATIVE=ON",
        f"-DCMAKE_PREFIX_PATH={os.environ.get('PREFIX', sys.prefix)}",
    ]
    # scikit-build reads these only from setup() arguments, not pyproject.toml.
    # Without them, files install under the environment root, not the package.
    options["packages"] = find_packages(include=["HiCOLA", "HiCOLA.*"])
    options["include_package_data"] = False
    options["package_data"] = {
        "HiCOLA": [
            path.relative_to("HiCOLA").as_posix()
            for path in Path("HiCOLA/Parameter_files").rglob("*")
            if path.suffix in {".ini", ".lua", ".txt", ".dat"}
        ],
    }

setup(version=version, **options)
