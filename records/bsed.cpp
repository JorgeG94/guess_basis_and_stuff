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
int main(int argc, char* argv[]) {
  //-------------------------------------//
  //-- change atom/basis set pair here --//
  //-------------------------------------//
  std::string atom(argv[1]);                                                    
  std::string basis(argv[2]);    

  //-- go into primary atom-basis set pair group --//
  H5::H5File file("bsed.h5", H5F_ACC_RDONLY);

  std::cout << "ATOM: " << atom << std::endl
            << "BASIS: " << basis << std::endl
            << std::endl;

  std::string atom_basis = atom + "/" + basis;
  H5::Group atom_basis_group = file.openGroup(atom_basis);

  //-- loop over shells in atom-basis set pair --//
  for (int shell = 1; shell <= atom_basis_group.getNumObjs(); ++shell) {
    //-- go into shell group --//
    std::string shell_str = std::to_string(shell);
    H5::Group shell_group = atom_basis_group.openGroup(shell_str);

    //-- create datasets for shell --//
    H5::DataSet shell_exponents = shell_group.openDataSet("Exponents");

    H5::DataSet shell_coefficients = shell_group.openDataSet("Coefficients");

    //-- define dataspaces for shell --//
    H5::DataSpace exp_dataspace = shell_exponents.getSpace();
    hsize_t num_exponents;
    int exp_rank = exp_dataspace.getSimpleExtentDims(&num_exponents, NULL);
    // std::cout << exp_rank << ": " << num_exponents << std::endl;

    H5::DataSpace coef_dataspace = shell_coefficients.getSpace();
    hsize_t num_coefficients; 
    int coef_rank = coef_dataspace.getSimpleExtentDims(&num_coefficients, NULL);
    // std::cout << coef_rank << ": " << num_coefficients[0] << ","
    //          << num_coefficients[1] << std::endl;

    //-- read from data set --//
    std::vector<double> shell_exponent_buf(num_exponents, 0.0);
    shell_exponents.read(shell_exponent_buf.data(),
                         H5::PredType::NATIVE_DOUBLE);

    std::vector<double> shell_coefficient_buf(num_coefficients, 0.0);
    shell_coefficients.read(shell_coefficient_buf.data(),
                         H5::PredType::NATIVE_DOUBLE);

    //-- print results --//
    std::cout << "SHELL " << shell
              << " EXPONENTS AND COEFFICIENTS:" << std::endl
              << "-----------------------------------------------------"
              << std::endl;
    for (int idx = 0; idx != num_exponents; ++idx) {
      if (num_coefficients != num_exponents) { // L shells
        std::cout << shell_exponent_buf[idx] << "     |     " <<
          shell_coefficient_buf[idx] << "     " <<
          shell_coefficient_buf[idx + num_exponents] << std::endl;
      } else {
        std::cout << shell_exponent_buf[idx] << "     |     "
                << shell_coefficient_buf[idx] << "     " << std::endl;
      }
    }
    std::cout << std::endl;
  }
  return 0;
}
