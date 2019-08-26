#include "silopp/silo++.hpp"
#include <iostream>

using namespace silo;

void show_dir(file_t& sf, std::string indent) {
  for ( auto arr : sf.toc_array() ) {
    std::cout << indent << arr << ", len = " << sf.var_length(arr) << std::endl;
  }
  for ( auto dir : sf.toc_dir() ) {
    std::cout << indent << "(D)" << dir << std::endl;
    sf.cd(dir);
    show_dir(sf, indent + "  " );
    sf.cd("..");
  }

}

int main(int argc, char** argv) {
  if ( argc != 2 ) {
    std::cout << "Incorrect number of arguments" << std::endl;
    return 0;
  }
  auto sf = open( std::string(argv[1]), Mode::Read );
  show_dir(sf,"");

  return 0;
}
