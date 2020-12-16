#include <array>
#include <cmath>
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
  std::string basis = "PCSeg-0";

  //-- go into primary atom-basis set pair group --//
  H5::H5File file("sadgss.h5", H5F_ACC_RDONLY);

  std::cout << "ATOM: " << atom << std::endl
            << "BASIS: " << basis << std::endl
            << std::endl;

  std::string atom_basis_pair = atom + "/" + basis;
  
  H5::DataSet guess_dataset = file.openDataSet(atom_basis_pair);

  H5::DataSpace guess_dataspace = guess_dataset.getSpace();
  hsize_t guess_num_elements;
  int guess_rank = guess_dataspace.getSimpleExtentDims(&guess_num_elements, NULL);
  int nbas = floor(std::sqrt(guess_num_elements));

  std::cout << guess_rank << ": " << guess_num_elements << std::endl;

  //-- read from data set --//
  std::vector<double> guess_buf(guess_num_elements, 0.0);
  guess_dataset.read(guess_buf.data(), H5::PredType::NATIVE_DOUBLE);

  //-- print results --//
  std::cout << "ATOM " << atom << " / BASIS SET " << basis
              << " SAD GUESS:" << std::endl
              << "-----------------------------------------------------"
              << std::endl 
              << "-----------------------------------------------------"
              << std::endl;
  std::cout << "ELEMENT     |     GUESS     " << std::endl; 
  std::cout << "-------     |     -----     " << std::endl; 
 
  for (int i = 0, idx = 0; i != nbas; ++i) {
    for (int j = 0; j != nbas; ++j, ++idx) {
      std::cout << i << "   " << j << "   |     " << guess_buf.at(idx) << std::endl; 
    }
  }   
  std::cout << std::endl;
  return 0;
}
