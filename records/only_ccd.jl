module BSEParse

using PyCall
using HDF5
using JSON

bse = pyimport("basis_set_exchange")
const am_to_shell_mapping = Dict(
    0 => "S", 
    1 => "P",
    2 => "D", 
    3 => "F", 
    4 => "G", 
    5 => "H", 
    6 => "I" 
  )

#===========================================================#
#= perform parsing operation on single atom/basis set pair =#
#===========================================================#
function parse_individual(atom::Dict{String,Any}, atomid::String, basis::String, 
  bsed::HDF5File)
  
  for (shell_num, shell) in enumerate(atom["electron_shells"])
    ang_mom = shell["angular_momentum"]
    shell_type_string = length(ang_mom) > 1 ? "L" : am_to_shell_mapping[ang_mom[1]]  
    
    exponents::Array{Float64,1} = parse.(Float64,shell["exponents"])
   
    coeff::Array{Float64} = begin
      if shell_type_string == "L"
        s_coeff = parse.(Float64,shell["coefficients"][1])
        p_coeff = parse.(Float64,shell["coefficients"][2])
        vcat(s_coeff,p_coeff)
      else
        parse.(Float64,shell["coefficients"][1])
      end
    end

    #println("Writing to $atomid/$basis/$shell_num...")
    h5write("bsed.h5",
      "$atomid/$basis/$shell_num/Shell Type", shell_type_string)
    h5write("bsed.h5",
      "$atomid/$basis/$shell_num/Exponents", exponents)
    h5write("bsed.h5",
      "$atomid/$basis/$shell_num/Coefficients", coeff)
    
    #=
    println("Exponents:")
    display(h5read("bsed.h5",
      "$atomid/$basis/$shell_num/Exponents"))
    println("")
    println("Contraction Coefficients:")
    display(h5read("bsed.h5",
      "$atomid/$basis/$shell_num/Coefficients"))
    println("")
    =#
  end
end

#===========================================================#
#= perform parsing operation on single atom/basis set pair =#
#===========================================================#
function parse_individual_cc_new(atom::Dict{String,Any}, atomid::String, basis::String, 
  bsed::HDF5File)

  if atomid == "Kr"
  
  # loop over the shells
  shell_num = 1
  for shell in atom["electron_shells"]
    ang_mom = shell["angular_momentum"]
    shell_type_string = am_to_shell_mapping[ang_mom[1]]

    println("#===============================================#")
    println("#== CREATE INITIAL EXPONENT/COEFFICIENT LISTS ==#")  
    println("#===============================================#")
    println()
    println("#== EXPONENTS ==#")
    exponents::Vector{Float64} = parse.(Float64,shell["exponents"])
    display(exponents); println()
    
    println()
    println("#== COEFFICIENTS ==#")
    coeffs::Vector{Vector{Float64}} = []
    coefs_01::Vector{Union{Nothing,Int64}} = []
    coefs_02::Vector{Union{Nothing,Int64}} = []
    for icoef in shell["coefficients"] 
      coeff::Vector{Float64} = parse.(Float64,icoef)
      display(coeff); println()
      push!(coeffs, coeff)
      push!(coefs_01, findfirst(i->(i == 1.0),coeff))
      push!(coefs_02, findfirst(i->(i == 0.0),coeff))
    end
    
    display(coefs_01); println()
    println(" ** ")
    display(coefs_02); println()
    
    indexcoefs_01 = findall(x -> x == nothing, coefs_01)
    indexcoefs_02 = findall(x -> x == length(coeffs[1]), coefs_02)
    if(length(indexcoefs_02) >  0)
      println("FOUND", indexcoefs_02)
    end
    #for i in indexcoefs_02 
    #  println("found at ", i)
    #end 


    println("#=====================================#")
    println("#== TRIM EXPONENT/COEFFICIENT LISTS ==#")  
    println("#=====================================#")
    idx_to_remove = findall(x->typeof(x)==Int64, coefs_01)
    idx_to_remove_d = findfirst(x->typeof(x)==Int64, coefs_02)
    println(" d to remove ", coefs_02[1])
    display(idx_to_remove)  
 
    println() 
    println("#== EXPONENTS ==#")
    exp_cc::Vector{Vector{Float64}} = [ deepcopy(exponents) ] 
    display(exp_cc); println()
    exps_to_remove::Vector{Vector{Float64}} = [] 
    for iexp_cc in exp_cc, idx in idx_to_remove
      println(iexp_cc)
      println("idx = ", idx, " / ", idx_to_remove)
      push!(exps_to_remove, [ iexp_cc[coefs_01[idx]]])
      #if (length(findall(x-> x == nothing, coefs_01)) == 0)
      if (length(indexcoefs_01) == 0)
        iexp_cc[coefs_01[idx]] = 0.0
        #println("no nothing")
      end
      if (length(indexcoefs_02) > 0 && idx == idx_to_remove[length(idx_to_remove)])
        iexp_cc[coefs_02[1]] = 0.0
      end
    end 
    for iexp in exps_to_remove
      push!(exp_cc, iexp)
    end
    display(exp_cc); println()
    for iexp_cc in exp_cc
      filter!(x->x!=0.0, iexp_cc)
    end 
    filter!(!isempty, exp_cc)
    display(exp_cc); println()
 
    println()
    println("#== COEFFICIENTS ==#")   
    coeff_cc::Vector{Vector{Float64}} = [ deepcopy(icoef) for icoef in coeffs ]
    coeffs_to_remove::Vector{Vector{Float64}} = [] 
    for (coef_idx, icoef_cc) in enumerate(coeff_cc), idx in idx_to_remove
      if typeof(coefs_01[coef_idx]) == Nothing 
        push!(coeffs_to_remove, [ icoef_cc[coefs_01[idx]] ])
        #icoef_cc[coefs_01[idx]] = 0.0
      end 
    end 
    for icoeff in coeffs_to_remove
      push!(coeff_cc, icoeff)
    end
    for icoeff_cc in coeff_cc
      filter!(x->x!=0.0, icoeff_cc)
    end 
    filter!(!isempty, coeff_cc)
    display(coeff_cc); println()
  
    exp_id = typeof(coefs_01[1]) == Int64 ? 0 : 1
    for coef_id in 1:length(coefs_01) 
      if typeof(coefs_01[coef_id]) == Int64
        exp_id += 1
      end
 
      h5write("bsed.h5",
        "$atomid/$basis/$shell_num/Shell Type", shell_type_string)
      h5write("bsed.h5",
        "$atomid/$basis/$shell_num/Exponents", exp_cc[exp_id])
      h5write("bsed.h5",
        "$atomid/$basis/$shell_num/Coefficients", coeff_cc[coef_id])
      
      shell_num += 1
    end #for icoef
  end #for shells
  
  
  println("#==========================================#")
  println("#== FINAL LIST OF EXPONENTS/COEFFICIENTS ==#") 
  println("#==========================================#")
  for ishl in 1:(shell_num-1)
    println("Shell Type $atomid/$basis/$ishl:")
    display(h5read("bsed.h5",
      "$atomid/$basis/$ishl/Shell Type"))
    println()
    println("Exponents $atomid/$basis/$ishl:")
    display(h5read("bsed.h5",
      "$atomid/$basis/$ishl/Exponents"))
    println()
    println("Contraction Coefficients $atomid/$basis/$ishl:")
    display(h5read("bsed.h5",
      "$atomid/$basis/$ishl/Coefficients"))
    println()
  end
  
  end # if
end #function

#===========================================#
#= parse all selected atom-basis set pairs =#
#===========================================#
function parse_all()
    #== open HDF5 file for writing ==#
    hdf5name::String = "bsed.h5"
    h5open(hdf5name, "w") do bsed
        flush(bsed)

        #== parse "sto" basis family ==#
        atoms::Array{String,1} = [
            "H", "He",
            "Li", "Be", "B", "C", "N", "O", "F", "Ne",
            "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
            "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
            "Ga", "Ge", "As", "Se", "Br", "Kr",
            "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
            "In", "Sn", "Sb", "Te", "I", "Xe",
            "Cs"
            ] #H-Cs

        #== parse correlation-consistent basis family ==#
        #basis_sets = ["cc-pVDZ", "aug-cc-pVDZ", "cc-pVTZ", "aug-cc-pVTZ"] 
        basis_sets = ["cc-pVTZ"] 
        
        for basis::String in basis_sets
            println("Basis: $basis")
            bs_dict = bse.get_basis(basis,fmt="json", header=false)
            bs_dict_json = JSON.parse(bs_dict)

            for atom in bs_dict_json["elements"] 
              parse_individual_cc_new(atom.second, atoms[parse(Int64,atom.first)], 
                basis, bsed)
            end
        end
        
    end
end

end

BSEParse.parse_all()
