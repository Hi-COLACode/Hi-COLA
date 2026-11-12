-- Shared settings for the small simulations.
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

local directory = debug.getinfo(1, "S").source:sub(2):match("^(.*[/\\])") or "./"
local repository = directory .. "../../"

-- Small periodic volume and deterministic initial conditions.
all_parameters_must_be_in_file = false
simulation_boxsize = 80.0
particle_Npart_1D = 16
particle_allocation_factor = 4.0
ic_nmesh = 16
force_nmesh = 32
pofk_nmesh = 32
ic_random_seed = 1234
ic_fix_amplitude = true
ic_initial_redshift = 49.0
ic_input_filename = repository .. "HiCOLA/Parameter_files/Backend_parameter_files/LCDM_1306dot3219_matterpower_z0.000.dat"

cosmology_model = "LCDM"
gravity_model = "GR"
cosmology_h = 0.7307
cosmology_Omegab = 0.02196 / cosmology_h^2
cosmology_OmegaCDM = 0.1274 / cosmology_h^2
cosmology_OmegaMNu = 0.0
force_linear_massive_neutrinos = false

output_particles = true
output_fileformat = "GADGET"
output_redshifts = {0.0}
timestep_nsteps = {5}
fof = false
pofk = true

-- Lightcone output options apply when a workflow enables lightcone generation.
lightcone = false
plc_pos_observer = {0.5, 0.5, 0.5}
plc_z_init = 0.25
plc_z_finish = 0.0
plc_use_replicas = false
plc_boundary_mpch = 0.0
plc_output_asciifile = true
plc_output_gadgetfile = false
plc_output_in_batches = false
plc_make_onion_density_maps = false
