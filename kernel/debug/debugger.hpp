#ifndef _DEBUGGER_HPP_
#define _DEBUGGER_HPP_
#include <string>
#include <vector>

#include "debug/nan.hpp"

namespace debug {
extern int timestep;
extern int world_rank;
extern int ens_label;
extern std::vector<double> dbls;
extern std::vector<int> ints;
extern std::vector<float> flts;
extern std::vector<std::string> strs;

inline void throw_error(std::string message) {
  throw std::runtime_error(message + " (rank " + std::to_string(world_rank) +
                           ")");
}
}  // namespace debug

#endif
