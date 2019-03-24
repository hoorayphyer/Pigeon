#ifndef  _PARTICLE_PROPERTIES_HPP_
#define  _PARTICLE_PROPERTIES_HPP_

#include "particle/map.hpp"

namespace particle {
  struct Properties {
    unsigned int mass_x = 0; // in terms of unit mass
    int charge_x = 0; // in terms of unit charge
    bool is_radiative = false;
  };
}

namespace particle {
  extern map<Properties> properties;
}

// namespace particle {
//   template <> constexpr unsigned int
//   mass_x<species::electron> = 1;
//   template <> constexpr unsigned int
//   mass_x<species::positron> = 1;

//   template <> constexpr unsigned int
//   charge_x<species::electron> = -1;
//   template <> constexpr unsigned int
//   charge_x<species::positron> = 1;
//   template <> constexpr unsigned int
//   charge_x<species::ion> = 1;
// }

// #include "./traits.hpp" // TODO this will break LAB
// namespace particle {
//   template <> constexpr unsigned int
//   mass_x<species::ion> = traits::ion_mass;
// }


#endif
