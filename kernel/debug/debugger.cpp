#include "debug/debugger.hpp"

namespace debug {
  int timestep {};
  int world_rank {};
  int ens_label = -1;
  std::vector<double> dbls {};
  std::vector<int> ints {};
  std::vector<float> flts {};
  std::vector<std::string> strs {};
}
