#ifndef  _PARTICLE_TRAITS_HPP_
#define  _PARTICLE_TRAITS_HPP_

namespace particle {
  enum class species : unsigned char
    { electron = 0, positron, ion, photon };

  enum class flag : unsigned int
    { empty = 0, secondary, ignore_force,
      ignore_deposit, annihilate, ignore_em,
      delimiter, traced };
}

namespace particle {
  template < species sp >
  constexpr unsigned int mass_x; // in terms of unit mass

  template < species sp >
  constexpr int charge_x; // in terms of unit charge

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
  mass_x<species::ion> = 5; // TODO move it somewhere else
  template <> constexpr unsigned int
  mass_x<species::photon> = 0;

  template <> constexpr unsigned int
  charge_x<species::electron> = -1;
  template <> constexpr unsigned int
  charge_x<species::positron> = 1;
  template <> constexpr unsigned int
  charge_x<species::ion> = 1;
  template <> constexpr unsigned int
  charge_x<species::photon> = 0;
}

#endif
