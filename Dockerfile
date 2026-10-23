FROM python:3.11-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends \
    g++ make git libfftw3-dev libfftw3-mpi-dev libgsl-dev liblua5.4-dev \
    openmpi-bin libopenmpi-dev \
    && rm -rf /var/lib/apt/lists/*
COPY requirements-test.txt /tmp/requirements-test.txt
RUN python -m pip install --no-cache-dir -r /tmp/requirements-test.txt
WORKDIR /workspace
CMD ["/bin/bash"]
