#ifndef _DEBUGGER_HPP_
#define _DEBUGGER_HPP_
#include "debug/nan.hpp"
#include <vector>
#include <string>

namespace debug {
  extern int timestep;
  extern std::vector<double> dbls;
  extern std::vector<int> ints;
  extern std::vector<float> flts;
  extern std::vector<std::string> strs;
}

#endif
