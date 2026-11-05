"""Exercise the installed frontend-to-COLASolver path on a small grid."""

import os
from pathlib import Path
import sys
import tempfile

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import sympy as sym

from HiCOLA import backend, USING_CPP
from HiCOLA.Frontend import expression_builder as eb
from HiCOLA.Utilities.Other import support

assert USING_CPP, "The native conda package must include COLASolver"
parameter_files = Path(backend.__file__).parent / "Parameter_files"
assert any(parameter_files.rglob("*.ini")), "The package must include parameter files"
Path(eb.__file__).resolve().relative_to(Path(sys.prefix).resolve())
assert sym.expand(eb.K_func(eb.X**3 * eb.phi**2, "X") - 3 * eb.X**2 * eb.phi**2) == 0

with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    # LCDM-limit tables exercise the same file interface as the frontend.
    a = np.geomspace(1e-4, 2.0, 512)
    e_squared = 0.3 / a**3 + 0.7
    support.write_data_flex(
        [a, np.sqrt(e_squared), -1.5 * 0.3 / a**3 / e_squared],
        root / "expansion.txt",
    )
    support.write_data_flex(
        [a, np.zeros_like(a), np.zeros_like(a), np.ones_like(a)],
        root / "force.txt",
    )
    k = np.geomspace(1e-5, 1e3, 512)
    np.savetxt(root / "power.txt", np.column_stack([k, 1e4 * k / (1 + (k / 0.1)**4)]))
    parameters = root / "run.lua"
    parameters.write_text('''
all_parameters_must_be_in_file = false
simulation_name = "CondaTest"
simulation_boxsize = 100.0
simulation_use_cola = false
cosmology_model = "HiCOLA"
cosmology_Omegab = 0.05
cosmology_OmegaCDM = 0.25
cosmology_OmegaK = 0.0
cosmology_Neffective = 3.046
cosmology_TCMB_kelvin = 2.7255
cosmology_h = 0.7
gravity_model = "HiCOLA"
gravity_model_screening = false
HiCOLA_expansion_filename = "expansion.txt"
HiCOLA_preforce_filename = "force.txt"
particle_Npart_1D = 8
ic_nmesh = 8
force_nmesh = 8
ic_initial_redshift = 9.0
ic_random_seed = 1234
ic_type_of_input = "powerspectrum"
ic_input_filename = "power.txt"
ic_input_redshift = 0.0
ic_sigma8_normalization = true
ic_sigma8 = 0.8
timestep_nsteps = {2}
output_folder = "."
output_redshifts = {0.0}
output_particles = false
pofk = true
pofk_nmesh = 8
''', encoding="utf-8")
    result = backend.run(parameters, cwd=root, capture_output=True, timeout=180)
    print(result.stdout)
    assert "Doing output 0 at redshift z = 0" in result.stdout
    spectra = list(root.glob("pofk*"))
    assert spectra, "COLASolver must write a power spectrum"
    for spectrum in spectra:
        values = np.loadtxt(spectrum)
        assert values.size and np.isfinite(values).all(), spectrum
