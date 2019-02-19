#ifndef  _TRAITS_HPP_
#define  _TRAITS_HPP_

#include "kernel/shape_predef.hpp"
#include "particle/pair_produce_predef.hpp"
#include "particle/species_predef.hpp"
#include "kernel/coordsys_predef.hpp"

namespace traits {
  using real_t = double;
  using ptc_state_t = unsigned long long;

  constexpr int DGrid = 2; // simulation dimension

  constexpr int DPtc = 3; // simulation dimension

  constexpr std::array<int, 3> grid_N { 10, 5, 2 };

  constexpr std::array<int, 3> grid_guard { 1, 1, 1 };

  constexpr std::array< std::array< real_t , 2 >, 3 >
  q_limit { 0.0, 100.0, 0.0, 50.0, 0.0, 20.0 }; // TODO initializer

  constexpr auto shape = knl::shape::Cloud_In_Cell;

  constexpr auto pair_produce_scheme = particle::PairScheme::Photon;

  constexpr auto coordinate_system = knl::coordsys::Cartesian;

  constexpr unsigned int ion_mass = 5;

  constexpr auto posion_inj = particle::species::positron; // posion = positron || ion in injection

  using deposit_j_t = long double;
};

#endif
