#== import modules here ==#
using PyCall
psi4 = pyimport("psi4")

#== psi4 setup ==# 
psi4.set_memory("500 MB")
psi4.core.set_output_file("w1.dat")

geom = psi4.geometry("
  O  -4.3997  1.0764  7.7009                                                   
  H  -3.9597  0.2664  7.3109                                                   
  H  -4.4297  0.9964  8.7009  
  symmetry=c1
")

options = Dict(
  "REFERENCE" => "RHF",
  "BASIS" => "PCSEG-0", 
  "SCF_TYPE"  => "DIRECT",
  "GUESS" => "SAD",
  "SAD_SCF_TYPE" => "DIRECT",
  "DF_SCF_GUESS" => false,
  "E_CONVERGENCE" => 1e-10,
  "D_CONVERGENCE" => 1e-10,
  "PRINT" => 4
)
psi4.set_options(options)

name = "scf/"*options["BASIS"]
scf_e, scf_wfn = psi4.energy(name, molecule=geom, return_wfn=true)

S = scf_wfn.S().to_array()
println("Overlap matrix:")
display(S); println()
