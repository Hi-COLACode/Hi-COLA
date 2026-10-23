-- Shared settings for the small simulations.
-- run_simulations.py appends particle/grid sizes, paths, model and lightcone settings.
-- Other optional parameters use ReadParameters.h defaults.

simulation_use_cola = true

-- Required background parameters and the reference primordial spectrum.
cosmology_OmegaK = 0.0
cosmology_Neffective = 3.046
cosmology_TCMB_kelvin = 2.7255
cosmology_As = 2.10521e-9
cosmology_ns = 0.9665

-- Use the reference random generator and a z=0 input matter spectrum.
ic_random_generator = "GSL"
ic_type_of_input = "powerspectrum"
ic_input_redshift = 0.0

-- Interlaced PCS power-spectrum estimator with shot-noise subtraction.
pofk_interlacing = true
pofk_subtract_shotnoise = true
pofk_density_assignment_method = "PCS"
