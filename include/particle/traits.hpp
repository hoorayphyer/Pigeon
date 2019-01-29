#ifndef  _PARTICLE_TRAITS_HPP_
#define  _PARTICLE_TRAITS_HPP_

namespace {
  enum class species : unsigned char
    { electron = 0, positron, ion, photon };

  enum class flag : unsigned int
    { empty = 0, secondary, ignore_force,
      ignore_deposit, annihilate, ignore_em,
      delimiter, traced };
}

namespace particle {
  template < species sp >
  constexpr unsigned int mass; // in terms of unit mass

  template < species sp >
  constexpr unsigned int charge; // in terms of unit charge

  template < species sp >
  constexpr bool is_radiative;
}

#endif
