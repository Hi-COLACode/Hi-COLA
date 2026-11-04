FROM python:3.11-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends \
    g++ make git libfftw3-dev libfftw3-mpi-dev libgsl-dev liblua5.4-dev \
    openmpi-bin libopenmpi-dev \
    && rm -rf /var/lib/apt/lists/*
COPY pyproject.toml /tmp/hicola/pyproject.toml
COPY HiCOLA /tmp/hicola/HiCOLA
RUN python -m pip install --no-cache-dir /tmp/hicola
WORKDIR /workspace
CMD ["/bin/bash"]
