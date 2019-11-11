#ifndef  _LOGGER_HPP_
#define  _LOGGER_HPP_

#include "logger/ofstream.hpp"
#include <iostream>
#include <sstream>

namespace lgr {
  extern ofstream<> file;
  extern ofstream<> debug;

  extern std::ostream& out;
  extern std::ostream& err;
  extern std::ostringstream str;
}

#endif
