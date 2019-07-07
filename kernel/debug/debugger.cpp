#include "debug/debugger.hpp"

namespace debug {
  int timestep {};
  int world_rank {};
  std::vector<double> dbls {};
  std::vector<int> ints {};
  std::vector<float> flts {};
  std::vector<std::string> strs {};
}
