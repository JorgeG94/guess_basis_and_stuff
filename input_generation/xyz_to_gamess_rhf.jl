include("gamess_custom_basis_sets/631+gd_split.jl")

const atom_to_atomic_number_mapping = Dict(
  "H " => 1.0,
  "He" => 2.0,
  "Li" => 3.0,
  "Be" => 4.0,
  "B " => 5.0,
  "C " => 6.0,
  "N " => 7.0,
  "O " => 8.0,
  "F " => 9.0,
  "Ne" => 10.0,
  "Na" => 11.0,
  "Mg" => 12.0,
  "Al" => 13.0,
  "Si" => 14.0,
  "P " => 15.0,
  "S " => 16.0,
  "Cl" => 17.0,
  "Ar" => 18.0,
  "K " => 19.0,
  "Ca" => 20.0,
  "Sc" => 21.0,
  "Ti" => 22.0,
  "V " => 23.0,
  "Cr" => 24.0,
  "Mn" => 25.0,
  "Fe" => 26.0,
  "Co" => 27.0,
  "Ni" => 28.0,
  "Cu" => 29.0,
  "Zn" => 30.0,
  "Ga" => 31.0,
  "Ge" => 32.0,
  "As" => 33.0,
  "Se" => 34.0,
  "Br" => 35.0,
  "Kr" => 36.0,
  "Rb" => 37.0,
  "Sr" => 38.0,
  "Y " => 39.0,
  "Zr" => 40.0,
  "Nb" => 41.0,
  "Mo" => 42.0,
  "Tc" => 43.0,
  "Ru" => 44.0,
  "Rh" => 45.0,
  "Pd" => 46.0,
  "Ag" => 47.0,
  "Cd" => 48.0,
  "In" => 49.0,
  "Sn" => 50.0,
  "Sb" => 51.0,
  "Te" => 52.0,
  "I " => 53.0,
  "Xe" => 54.0
)

function create_input(input)
  xyz_array::Vector{String} = []
  xyz = open(input, "r") do file
    fragment = readlines(file)  
    xyz_array = fragment[3:end]
  end
  natoms = length(xyz_array)

  output = open(joinpath(@__DIR__, input[1:(end-4)]*".inp"), "w") do file
    write(file, " \$scf dirscf=.t. damp=.t. shift=.f. \$end\n")
    write(file, "   soscf=.f. diis=.t. fdiff=.f. extrapolate=.f. \$end\n")
    write(file, " \$system mwords=500 memddi=500 \$end\n")
    write(file, " \$intgrl schwarz=.t. intomp=1 shfock=.f. \$end\n")
    write(file, " \$guess guess=huckel \$end\n")
    write(file, " \$contrl scftyp=rhf runtyp=energy nprint=-5\n")
    write(file, "     icharg=0 ispher=0 inttyp=eric \n")
    write(file, "     maxit=50 \$end\n")
    write(file, " \$end\n")
    write(file, " \$data\n")
    for line in xyz_array 
      atom = line[1:2]
      newline = line[1:2]*"   $(atom_to_atomic_number_mapping[atom])   "*
        line[3:end]
      write(file, "$newline\n")
      
      newline = six_three_one_plus_g_d_split[atom]
      write(file, "$newline\n")
    end
    write(file, " \$end\n")
  end 
end

create_input(ARGS[1])
