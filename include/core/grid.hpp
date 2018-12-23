#ifndef  _GRID_HPP_
#define  _GRID_HPP_

#include "traits.hpp"

template < int I_Dim >
struct grid {
private:
  static constexpr Real
  lower = I_Dim < traits::Dgrid ? traits::q_limit[I_Dim][0] : 0;

  static constexpr Real
  upper = I_Dim < traits::Dgrid ? traits::q_limit[I_Dim][1] : 1;

public:
  static constexpr int
  guard = I_Dim < traits::Dgrid ? traits::grid_guard[I_Dim] : 0;

  static constexpr Real
  delta = I_Dim < traits::Dgrid ?
                  ( traits::q_limit[I_Dim][1] - traits::q_limit[I_Dim][0] ) / traits::grid_N[I_Dim]
                  : 1;

  // IsCentered = true corresponds to unstaggered
  template < bool IsCentered >
  static constexpr
  Real abscissa( const int i ) noexcept {
    if constexpr ( I_Dim < traits::Dgrid )
      return lower - delta * guard + ( IsCentered ? 0.5 : 1.0 ) * delta + i * delta;
    else
      return lower;
  }

};


template < int I_Dim >
struct grid_patch {
  static constexpr int
  extent = I_Dim < traits::Dgrid ? traits::grid_N[I_Dim] / traits::partition[I_Dim] : 1; // not including guard cells

  static int origin; // the equivalent of local index 0 in the global grid

  static constexpr int extent_with_guard() {
    return extent + 2 * grid<I_Dim>::guard;
  }
};

#endif
