local directory = debug.getinfo(1, "S").source:sub(2):match("^(.*[/\\])") or "./"
dofile(directory .. "backend_minimal.lua")

simulation_name = "gr"
output_folder = directory .. "../../Hi-COLA_Output/gr/output"
