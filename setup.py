"""Build the Python frontend, optionally including the COLASolver executable."""

import os
import re
import sys

native = os.environ.get("HICOLA_BUILD_NATIVE", "0") == "1"
if native:
    from skbuild import setup
else:
    from setuptools import setup

version = os.environ.get("PKG_VERSION", "0.0.0")
if os.environ.get("GITHUB_REF_TYPE") == "tag":
    tag = os.environ.get("GITHUB_REF_NAME", "")
    if not re.fullmatch(r"v\d+\.\d+\.\d+", tag):
        raise ValueError("Release tags must have the form vMAJOR.MINOR.PATCH")
    version = tag[1:]

options = {}
if native:
    options["cmake_args"] = [
        "-DHICOLA_BUILD_NATIVE=ON",
        f"-DCMAKE_PREFIX_PATH={os.environ.get('PREFIX', sys.prefix)}",
    ]

setup(version=version, **options)
