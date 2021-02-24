#== import modules here ==#
using PyCall
bse = pyimport("basis_set_exchange")

using HDF5

#== initialization ==#
upper_to_lower = Dict(
                      "A" => "a",
                      "B" => "b",
                      "C" => "c",
                      "D" => "d",
                      "E" => "e",
                      "F" => "f",
                      "G" => "g",
                      "H" => "h",
                      "I" => "i",
                      "J" => "j",
                      "K" => "k",
                      "L" => "l",
                      "M" => "m",
                      "N" => "n",
                      "O" => "o",
                      "P" => "p",
                      "Q" => "q",
                      "R" => "r",
                      "S" => "s",
                      "T" => "t",
                      "U" => "u",
                      "V" => "v",
                      "W" => "w",
                      "X" => "x",
                      "Y" => "y",
                      "Z" => "z"
                     )
#==============================#
#== get atom/basis set pairs ==# 
#==============================#

basis_list = Vector{String}([])
status_bsed = h5open("../records/bsed.h5","r") do bsed
  #== loop over atoms ==#
  for symbol in names(bsed)
    #== loop over basis sets in atom ==#
    for basis in names(bsed["$symbol"])
      #== add basis set to pair list if unique ==#
      if !in(basis, basis_list) 
        push!(basis_list, basis)
      end
    end
  end
end

for basis in basis_list
  
  bs_psi4 = bse.get_basis(basis,fmt="psi4", header=false)
  bs_psi4 = replace(bs_psi4, "D+" => "E+")       
  bs_psi4 = replace(bs_psi4, "D-" => "E-")       
 
  bs_name = basis
  bs_name = replace(bs_name, "*" => "s")       
  bs_name = replace(bs_name, "+" => "p")       
  bs_name = replace(bs_name, "(" => "_")       
  bs_name = replace(bs_name, ")" => "_")       
  bs_name = replace(bs_name, "," => "_")  
  for letter in upper_to_lower
    bs_name = replace(bs_name, letter[1] => letter[2])  
  end 

  println(bs_name)    

  bs_file = open("$bs_name.gbs", "w") do file
    write(file, bs_psi4)    
  end
end

