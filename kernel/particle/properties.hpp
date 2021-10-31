#pragma once

#include <string>

namespace particle {
struct Properties {
  float mass_x = 0;    // in terms of unit mass
  float charge_x = 0;  // in terms of unit charge
  std::string name = "";
  std::string nickname = "";
};
}  // namespace particle
