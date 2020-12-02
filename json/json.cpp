#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main() {
  //-- create json object --//
  std::ifstream ifs("fig1a.json");
  std::string content( (std::istreambuf_iterator<char>(ifs) ),
                       (std::istreambuf_iterator<char>()    ) );

  nlohmann::json input_file = nlohmann::json::parse(content);
  //std::cout << input_file.dump(2) << std::endl << std::endl;

  //-- parse molecular structure --//
  auto molecule = input_file.at("molecule");

  auto coords = molecule.at("geometry").get<std::vector<double> >();
  auto symbols = molecule.at("symbols").get<std::vector<std::string> >();
  
  std::cout << "------------------------------------" << std::endl;
  std::cout << "   Parsing molecular structure...   " << std::endl;
  std::cout << "------------------------------------" << std::endl;
  std::cout << "Atom  |         x         y         z        " << std::endl;
  for (int isymbol = 0; isymbol != symbols.size(); ++isymbol) { 
    std::cout << symbols[isymbol] << "    =>    " << coords[3*isymbol] << ", " 
      << coords[3*isymbol + 1] << ", " << coords[3*isymbol + 2] << std::endl;
  }
  std::cout << std::endl;
  
  //-- parse calculation information --//
  auto model = input_file.at("model");

  auto basis = model.at("basis").get<std::string>(); 
  auto method = model.at("method").get<std::string>();

  std::cout << "----------------------------------------" << std::endl;
  std::cout << "   Parsing calculation information...   " << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "The calculation type we are running is: " << method << std::endl;
  std::cout << "The basis set we are using is: " << basis << std::endl; 
  std::cout << std::endl;

  // -- parse keywords --// 
  auto keywords = input_file.at("keywords");
  
  auto scf = keywords.at("scf");
  
  std::cout << "-------------------------------------" << std::endl;
  std::cout << "   Parsing calculation keywords...   " << std::endl;
  std::cout << "-------------------------------------" << std::endl;
  
  std::cout << "We have set the following SCF keywords: " << std::endl;
  for (nlohmann::json::iterator it = scf.begin(); it != scf.end(); ++it) { 
    std::cout << it.key() << " => " << it.value() << std::endl; 
  } 
  std::cout << std::endl;
  
  //-- we are done! --//
  return 0;
}

