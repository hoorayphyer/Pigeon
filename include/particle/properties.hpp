#ifndef  _PARTICLE_PROPERTIES_HPP_
#define  _PARTICLE_PROPERTIES_HPP_

#include "species_predef.hpp"

namespace particle {
  template < species sp >
  constexpr unsigned int mass_x = 0; // in terms of unit mass

  template < species sp >
  constexpr int charge_x = 0; // in terms of unit charge

  template < species sp >
  constexpr bool is_radiative = ( sp == species::electron || sp == species::positron );

  template < species sp >
  constexpr bool is_charged = (charge_x<sp> != 0);

  template < species sp >
  constexpr bool is_massive = (mass_x<sp> != 0);
}

namespace particle {
  template <> constexpr unsigned int
  mass_x<species::electron> = 1;
  template <> constexpr unsigned int
  mass_x<species::positron> = 1;

  template <> constexpr unsigned int
  charge_x<species::electron> = -1;
  template <> constexpr unsigned int
  charge_x<species::positron> = 1;
  template <> constexpr unsigned int
  charge_x<species::ion> = 1;
}

#include "traits.hpp"
namespace particle {
  template <> constexpr unsigned int
  mass_x<species::ion> = traits::ion_mass;
}


#endif
