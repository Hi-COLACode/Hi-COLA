local directory = debug.getinfo(1, "S").source:sub(2):match("^(.*[/\\])") or "./"
dofile(directory .. "backend_minimal.lua")

simulation_name = "cubic_galileon"
output_folder = directory .. "../../Hi-COLA_Output/cubic_galileon/output"

cosmology_model = "HiCOLA"
gravity_model = "HiCOLA"
HiCOLA_expansion_filename = directory .. "../../Hi-COLA_Output/frontend/cubic_galileon/cubic_galileon_WMAP_expansion.txt"
HiCOLA_preforce_filename = directory .. "../../Hi-COLA_Output/frontend/cubic_galileon/cubic_galileon_WMAP_force.txt"
gravity_model_screening = true
gravity_model_HiCOLA_smoothing_filter = "gaussian"
gravity_model_HiCOLA_smoothing_scale_over_boxsize = 3.0 / simulation_boxsize
gravity_model_screening_enforce_largescale_linear = true
gravity_model_screening_linear_scale_hmpc = 0.1
