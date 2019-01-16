#ifndef  _TRAITS_HPP_
#define  _TRAITS_HPP_

#include "types.hpp"


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


};


constexpr Species Ion = { "Ion", 1, 5, false };

#endif
