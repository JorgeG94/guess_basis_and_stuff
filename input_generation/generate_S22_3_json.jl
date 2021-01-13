include("xyz_to_json.jl")

directory = joinpath(@__DIR__) * "/S22_3/"                                          
all_files = readdir(directory)                                                
all_files .= directory .* all_files                                           
#display(all_files)

create_input_rhf.(all_files)
