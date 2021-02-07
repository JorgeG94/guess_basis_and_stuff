include("xyz_to_gamess_rhf.jl")

directory = joinpath(@__DIR__) * "/S22_3/"                                          
all_files = readdir(directory)                                                
all_files .= directory .* all_files                                           
all_files = filter(x->occursin(".xyz",x), all_files)

create_input.(all_files)
