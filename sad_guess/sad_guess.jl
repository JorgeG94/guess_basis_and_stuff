#== import modules here ==#
using PyCall
psi4 = pyimport("psi4")

using HDF5

#== psi4 setup ==# 
psi4.set_memory("500 MB")
psi4.core.set_output_file("sad_guess.dat")

#==============================#
#== get atom/basis set pairs ==# 
#==============================#
pairs = Vector{String}([])
status_bsed = h5open("bsed.h5","r") do bsed
  #== loop over atoms ==#
  for symbol in names(bsed)
    #== loop over basis sets in atom ==#
    for basis in names(bsed["$symbol"])
      #== add atom/basis-set pair to pair list ==#
      pair = "$symbol/$basis"
      push!(pairs, pair)
    end
  end
end

status_sadgss = h5open("sadgss.h5","w") do sadgss
  flush(sadgss)
  for pair in pairs
    pair_regex = match(r"(.*)/(.*)", pair)
    symbol = pair_regex.captures[1]
    basis = pair_regex.captures[2]
    
    println("$symbol/$basis")
    
    try 
      geom = psi4.geometry("""
        $symbol  0.0 0.0 0.0         
      """)
      
      options = Dict(
        "REFERENCE" => "UHF",
        "BASIS" => basis, 
        "SCF_TYPE"  => "PK",
        "E_CONVERGENCE" => 1e-10,
        "D_CONVERGENCE" => 1e-10
      )
      psi4.set_options(options)

      name = "scf/"*options["BASIS"]
      scf_e, scf_wfn = psi4.energy(name, molecule=geom, return_wfn=true)

      density_a = scf_wfn.Da().array_interface()[1]
      #density_b = scf_wfn.Db().array_interface()[1]
      display(density_a[:])
       
      write(sadgss, pair, density_a[:])
    catch
      continue
    end
  end
end
