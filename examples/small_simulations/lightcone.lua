local directory = debug.getinfo(1, "S").source:sub(2):match("^(.*[/\\])") or "./"
dofile(directory .. "backend_minimal.lua")

simulation_name = "lightcone"
output_folder = directory .. "../../Hi-COLA_Output/lightcone/output"

-- Replicate the box to cover the lightcone radial extent.
lightcone = true
plc_use_replicas = true
