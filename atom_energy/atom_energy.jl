#== import modules here ==#
using PyCall
psi4 = pyimport("psi4")

using HDF5
#using MPI

#== psi4 setup ==# 
psi4.set_memory("500 MB")
psi4.core.set_output_file("atom_energy.dat")

#===============#
#== MPI setup ==#
#===============#
#MPI.Init()
#comm = MPI.COMM_WORLD

#mpi_rank = MPI.Comm_rank(comm)
#mpi_size = MPI.Comm_size(comm)

#==============================#
#== get atom/basis set pairs ==# 
#==============================#

pairs = Vector{String}([])
status_bsed = h5open("../records/bsed.h5","r") do bsed
  #== loop over atoms ==#
  for symbol in keys(bsed)
    #== loop over basis sets in atom ==#
    for basis in keys(bsed["$symbol"])
      #== add atom/basis-set pair to pair list ==#
      pair = "$symbol/$basis"
      push!(pairs, pair)
    end
  end
end

status_eatom = h5open("eatom.h5","w") do eatom
  #if mpi_rank == 0
    flush(eatom)
  #end
  
  #symbol = "O" 
  #basis = "6-31G" 
  #pair = "$symbol/$basis"
  
  #MPI.Barrier(comm)
  #== mp2 energies ==#
  for (pair_idx, pair) in enumerate(pairs)
    #if mpi_rank != (pair_idx-1)%mpi_size continue end

    pair_regex = match(r"(.*)/(.*)", pair)
    symbol = pair_regex.captures[1]
    basis = pair_regex.captures[2]

    if symbol == "H"  
    println("$symbol/$basis")
    
    try 
      geom = psi4.geometry("""
        $symbol  0.0 0.0 0.0         
        symmetry=c1
      """)
      
      options = Dict(
        "REFERENCE" => "ROHF",
        "BASIS" => basis, 
        "SCF_TYPE"  => "PK",
        "GUESS" => "SAD",
	      "PUREAM" => false,
        #"SAD_SCF_TYPE" => "DIRECT",
        "E_CONVERGENCE" => 1e-10,
        "D_CONVERGENCE" => 1e-10,
        #"MAXITER" => 1,
        "FAIL_ON_MAXITER" => false,
        "MP2_TYPE" => "CONV",
        "QC_MODULE" => "OCC"
      )
      psi4.set_options(options)

      mp2_name = "mp2/"*options["BASIS"]
      mp2_e, mp2_wfn = psi4.energy(mp2_name, molecule=geom, return_wfn=true)
      write(eatom, "RIMP2/"*pair, mp2_e)
    catch e                                                                       
      bt = catch_backtrace()                                                      
      msg = sprint(showerror, e, bt)                                              
      println(msg)       
      continue
    end
    end
  end

  #== rhf energies ==#
  for (pair_idx, pair) in enumerate(pairs)
    #if mpi_rank != (pair_idx-1)%mpi_size continue end

    pair_regex = match(r"(.*)/(.*)", pair)
    symbol = pair_regex.captures[1]
    basis = pair_regex.captures[2]
 
    println("$symbol/$basis")
    
    if symbol == "H"
    try 
      geom = psi4.geometry("""
        $symbol  0.0 0.0 0.0         
        symmetry=c1
      """)
      
      options = Dict(
        "REFERENCE" => "ROHF",
        "BASIS" => basis, 
        "SCF_TYPE"  => "PK",
        "GUESS" => "SAD",
	      "PUREAM" => false,
        #"SAD_SCF_TYPE" => "DIRECT",
        "E_CONVERGENCE" => 1e-10,
        "D_CONVERGENCE" => 1e-10,
        #"MAXITER" => 1,
        "FAIL_ON_MAXITER" => false
      )
      psi4.set_options(options)

      scf_name = "scf/"*options["BASIS"]
      scf_e, scf_wfn = psi4.energy(scf_name, molecule=geom, return_wfn=true)
      write(eatom, "RHF/"*pair, scf_e)
    catch e                                                                       
      bt = catch_backtrace()                                                      
      msg = sprint(showerror, e, bt)                                              
      println(msg)       
      continue
    end
    end
  end
end

#MPI.Finalize()
