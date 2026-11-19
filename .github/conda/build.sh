#!/bin/bash
# Build and install the serial COLASolver into the environment, then the Python frontend.
set -euxo pipefail
cmake -S . -B build -G Ninja ${CMAKE_ARGS} \
    -DCMAKE_BUILD_TYPE=Release \
    -DHICOLA_USE_HEALPIX=ON \
    -DCMAKE_PREFIX_PATH="${PREFIX}" \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}"
cmake --build build --parallel "${CPU_COUNT}"
cmake --install build
"${PYTHON}" -m pip install . --no-deps --no-build-isolation -vv
