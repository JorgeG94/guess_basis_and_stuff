Base.include(@__MODULE__, "../../get_indat.jl")
#=
function xyz_to_geometry(xyzfile)
  input_file = []
  #println(joinpath(@__DIR__, xyzfile))
  open(joinpath(@__DIR__, xyzfile)) do file
    input_file = readlines(file)
  end

  natoms = parse(Int, input_file[1])
  #println(natoms)
  geometry = input_file[3:end]
  #display(geometry); println()

  coords_vector = Vector{Float64}([])
  symbols = Vector{String}([])
  for geometry_line in geometry
    coords_match = eachmatch(r"([-]?[\d]{1,}\.[\d]{1,})",geometry_line)
    atom_coords = collect(coords_match)
    atom_coords_match = map(x -> x.match, atom_coords) 
    append!(coords_vector, parse.(Float64,atom_coords_match[1:3]))

    symbol_match = match(r"[A-Za-z]{1,2}",geometry_line)
    push!(symbols, symbol_match.match)
  end
  #coords = reshape(coords_vector, (3,length(geometry)))
  #coords = transpose(coords)
  coords = coords_vector

  #@assert natoms == size(coords)[1]
  #@assert natoms == length(symbols)

  #display(coords); println()
  #display(symbols); println()

  charge = 0

  return coords, symbols, charge
end
=#
function create_input_mbe()
  #== read in fragment files ==#
  directory = joinpath(@__DIR__) * "/"
  all_files = readdir(directory)
  all_files .= directory .* all_files 
  
  #input_match = match.(r"frag[0-9]+.xyz",all_files)
  #input_match_filter = filter!(x -> x != nothing, input_match)
  #inputs = map(x->x.match, input_match)

  num_frags = begin
    nfrags = 0
    Nfrag = open("../Nfrags", "r") do file
      argv = readlines(file)
      nfrags = parse(Int64, argv[2])
    end 
    nfrags
  end
 
  fragments = [ "$i" for i in 1:num_frags ]
  inputs = "frag" .* fragments .* ".xyz"
  display(inputs)
  
  fmobnd::Vector{String} = []
  for frag in 1:num_frags
    bnd = open("../BROKENBONDS", "r") do file
      fmobnd = readlines(file)  
      #push!(fmobnd, bonds[2:end])
    end
  end
 
  coords, symbols, indat = get_indat("../molecule.xyz", num_frags)
  charge = 0
 
  #== write input json file ==#
  output = open("input.json", "w") do file
    write(file, "{\n")
    #== write in molecule information ==#
    write(file, "  \"molecule\": {\n")
    #== write fragment geometry ==#
    write(file, "    \"geometry\" : [\n") 
    for icoord in 1:length(coords)
      write(file, "      ") 
          
      value = coords[icoord]
      if icoord == length(coords)
        write(file, "$value\n")
      elseif icoord%3 == 0 
        write(file, "$value,\n")
      else
        write(file, "$value,")
      end
    end
    write(file, "    ],\n") 
   
    #== write fragment symbols ==#
    write(file, "    \"symbols\" : [\n") 
    for iatom in 1:length(symbols)
      write(file, "      ") 
      symbol = symbols[iatom]
      if iatom == length(symbols)
        write(file, "\"$symbol\"\n")
      elseif iatom%5 == 0 
        write(file, "\"$symbol\",\n")
      else
        write(file, "\"$symbol\",")
      end
    end
    write(file, "    ],\n")
 
    #== assign atoms to fragments ==#
    write(file, "    \"fragments\" : {\n") 
    write(file, "      \"nfrag\": $num_frags,\n")
    write(file, "      \"fragid\": [\n") 
    for ifrag in 1:length(indat)
      write(file, "      ") 
      fragment = indat[ifrag]
      if ifrag == length(indat) 
        write(file, "$fragment\n")
      elseif ifrag%5 == 0
        write(file, "$fragment,\n")
      else
        write(file, "$fragment,")
      end
    end
    write(file, "      ],\n")

    #== write broken bonds ==#
    write(file, "      \"broken_bonds\" : [\n") 
    for ibnd in fmobnd[2:end]
      bond_match = eachmatch(r"[0-9]+",ibnd)
      broken_bonds = collect(bond_match)
      broken_bonds_match = map(x -> x.match, broken_bonds)
      if ibnd == fmobnd[end] 
        write(file,"        $(broken_bonds_match[1]), $(broken_bonds_match[2])\n")
      else
        write(file,"        $(broken_bonds_match[1]), $(broken_bonds_match[2]),\n")
      end
    end
    write(file, "      ],\n")
     
    #== write fragment charge ==#
    #write(file, "    \"molecular_charge\": $charge\n") 
    #write(file, "  },\n")
    
    #== write fragment charges ==#
    write(file, "      \"fragment_charges\" : [\n") 
    for ifrag in 1:num_frags
      write(file, "        ") 
      fragment = indat[ifrag]
      if ifrag == num_frags 
        write(file, "0\n")
      elseif ifrag%5 == 0
        write(file, "0,\n")
      else
        write(file, "0,")
      end
    end
    write(file, "      ]\n")
    write(file, "    }\n")
    write(file, "  },\n")
  
    #== write calculation driver and model information ==# 
    write(file,"  \"driver\": \"energy\",\n")
    write(file,"  \"model\": {\n")
    write(file,"    \"method\": \"MBE-RIMP2\",\n")
    write(file,"    \"basis\": \"6-31G\",\n")
    write(file,"    \"aux_basis\": \"cc-pVDZ\"\n")
    write(file,"  },\n")
    
    #== write keywords ==# 
    write(file,"  \"keywords\": {\n")
    write(file,"    \"scf\": {\n")
    write(file,"      \"niter\": 200\n")
    write(file,"    },\n")
    write(file,"    \"rimp2\": {\n")
    write(file,"      \"box_dim\": 15\n")
    write(file,"    }\n")
    write(file,"  }\n")
    write(file,"}")
  end 
end

function create_input_rhf(input)
  #== write input json file ==#
  filename = input[1:(end-4)]
  #display(filename)
  output = open("$filename.json", "w") do file
    write(file, "{\n")
    #== write in molecule information ==#
    write(file, "  \"molecule\": {\n")
    #== write each fragment's molecule information ==#

    coords, symbols, charge = xyz_to_geometry(input)
    #== write fragment geometry ==#
    write(file, "    \"geometry\" : [\n") 
    for icoord in 1:length(coords)
      write(file, "      ") 
        
      value = coords[icoord]
      #println(value)
      if icoord == length(coords)
        write(file, "$value\n")
      elseif icoord%3 == 0 
        write(file, "$value,\n")
      else
        write(file, "$value,")
      end
    end
    write(file, "    ],\n") 
 
    #== write fragment symbols ==#
    write(file, "    \"symbols\" : [\n") 
    for iatom in 1:length(symbols)
      write(file, "      ") 
      symbol = symbols[iatom]
      if iatom == length(symbols)
        write(file, "\"$symbol\"\n")
      elseif iatom%5 == 0 
        write(file, "\"$symbol\",\n")
      else
        write(file, "\"$symbol\",")
      end
    end
    write(file, "    ],\n") 
    #== write fragment charge ==#
    write(file, "    \"molecular_charge\": $charge\n") 
    write(file, "  },\n")
 
    #== write calculation driver and model information ==# 
    write(file,"  \"driver\": \"energy\",\n")
    write(file,"  \"model\": {\n")
    write(file,"    \"method\": \"RHF\",\n")
    write(file,"    \"basis\": \"6-31+G*\"\n")
    write(file,"  },\n")
    
    #== write keywords ==# 
    write(file,"  \"keywords\": {\n")
    write(file,"    \"scf\": {\n")
    write(file,"      \"guess\":\"hcore\",\n")
    write(file,"      \"dele\":1.0, \n")
    write(file,"      \"rmsd\":2.0E-5 \n")
    write(file,"    }\n")
    write(file,"  }\n")
    write(file,"}")
  end 
end

create_input_mbe()
