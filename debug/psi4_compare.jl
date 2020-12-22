#== import modules here ==#
using PyCall
psi4 = pyimport("psi4")

using HDF5

#== psi4 setup ==# 
psi4.set_memory("500 MB")
psi4.core.set_output_file("psi4_compare.dat")

geom = psi4.geometry("
  O      1.00339        0.03093        0.09503                                
  O      2.32339        0.03093        0.09503                                
  H      0.68006        -0.87208        -0.04962                              
  H      2.64672        0.93395        0.23969                                 
  symmetry=c1
")

options = Dict(
  "REFERENCE" => "RHF",
  "BASIS" => "PCSeg-1", 
  "SCF_TYPE"  => "DIRECT",
  "GUESS" => "SAD",
  "SAD_SCF_TYPE" => "DIRECT",
  "DF_SCF_GUESS" => false,
  "PUREAM" => false,
  "E_CONVERGENCE" => 1e-10,
  "D_CONVERGENCE" => 1e-10,
  "PRINT" => 4
)
psi4.set_options(options)

name = "scf/"*options["BASIS"]
scf_e, scf_wfn = psi4.energy(name, molecule=geom, return_wfn=true)

S_psi4 = scf_wfn.S().to_array()

#== standalone setup ==#
debug = h5open("debug.h5", "r")

tmp = debug["RHF"]["Iteration-None"]["S"][:]
S_standalone = reshape(tmp,(Int64(sqrt(length(tmp))), Int64(sqrt(length(tmp)))))
#display(S_standalone)

close(debug)

#== comparison ==#
S_diff = S_standalone .- S_psi4
#display(S_standalone[:])
for i in 1:length(S_diff)
  if S_diff[i] > 0.00001
    println(i," ", S_standalone[i]," ",  S_psi4[i]," ",  S_diff[i])
  end
end

leng = Int64(sqrt(length(S_diff)))

for i in 0:(leng-1)
  index = (leng+1)*i + 1
  println(index, " ", S_standalone[index], " ", S_psi4[index], " ", S_diff[index])
end
