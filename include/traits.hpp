#ifndef  _TRAITS_HPP_
#define  _TRAITS_HPP_

#include "types.hpp"
#include <experimental/array>

namespace std {
  using experimental::make_array;
}


struct traits {
  static constexpr int
  Dgrid = 2; // simulation dimension

  static constexpr std::array<int, 3>
  grid_N { 10, 5, 2 };

  static constexpr std::array<int, 3>
  grid_guard { 1, 1, 1 };

  static constexpr std::array< std::array< Real , 2 >, 3 >
  q_limit { 0.0, 100.0, 0.0, 50.0, 0.0, 20.0 }; // TODO initializer

  static constexpr ct_string
  shape {"Cloud In Cell"};

  static constexpr auto
  active_species = std::make_array
    ( "electron", "positron", "photon" );

};

// TODO ion charge set elsewhere

#endif
