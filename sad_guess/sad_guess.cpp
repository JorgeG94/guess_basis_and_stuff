#include <array>
#include <iostream>
#include <string>
#include <vector>

#include <H5Cpp.h>

//------------------------------------------------------------------//
// This code loops over the shells of a specified atom with a given //
// basis set and prints out the exponents for each shell, acquiring //
// basis set information from JuliaChem.jl's bsed.h5 basis set      //
// database. Coefficients will be arriving soon.                    //
//------------------------------------------------------------------//
int main() {
  //-------------------------------------//
  //-- change atom/basis set pair here --//
  //-------------------------------------//
  std::string atom = "O";
  std::string basis = "6-31G";

  //-- go into primary atom-basis set pair group --//
  H5::H5File file("sadgss.h5", H5F_ACC_RDONLY);

  std::cout << "ATOM: " << atom << std::endl
            << "BASIS: " << basis << std::endl
            << std::endl;

  std::string atom_basis_pair = atom + "/" + basis;
  H5::DataSet pair_dataset = file.openDataSet(atom_basis_pair);

  H5::DataSpace pair_dataspace = pair_dataset.getSpace();
  hsize_t num_elements;
  int pair_rank = pair_dataspace.getSimpleExtentDims(&num_elements, NULL);
  std::cout << pair_rank << ": " << num_elements << std::endl;

  //-- read from data set --//
  std::vector<double> density_buf(num_elements, 0.0);
  pair_dataset.read(density_buf.data(), H5::PredType::NATIVE_DOUBLE);

  //-- print results --//
  std::cout << "ATOM " << atom << " / BASIS SET " << basis
              << " PAIR DENSITY:" << std::endl
              << "-----------------------------------------------------"
              << std::endl;
  for (int idx = 0; idx != num_elements; ++idx) {
    std::cout << density_buf[idx] << std::endl; 
  }
  std::cout << std::endl;
  return 0;
}
