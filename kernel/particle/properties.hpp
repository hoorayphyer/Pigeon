#ifndef _PARTICLE_PROPERTIES_HPP_
#define _PARTICLE_PROPERTIES_HPP_

#include "particle/map.hpp"
#include <string>

namespace particle {
  struct Properties {
    unsigned int mass_x = 0; // in terms of unit mass
    int charge_x = 0; // in terms of unit charge
    std::string name = "";
  };

  // defined at particle updater
  extern map<Properties> properties;
}

#endif
