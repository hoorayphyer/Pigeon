#ifndef  _LOGGER_HPP_
#define  _LOGGER_HPP_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

namespace lgr {
  extern std::string indent;

  extern std::fstream file;
  extern std::fstream debug;

  auto& out = std::cout;
  auto& err = std::cerr;
  extern std::stringstream str;
}

#endif
