#ifndef  _GRID_HPP_
#define  _GRID_HPP_

#include "types.hpp"
#include <array>

template < typename T >
struct Gridline {
  T lower;
  T upper;
  int guard;
  T delta;

  static_assert( std::is_floating_point_v<T> );

  constexpr Gridline() noexcept: Gridline( 0.0, 1.0, 0, 1 ) {}

  constexpr Gridline( T grid_lower, T grid_upper, int grid_guard, int grid_N ) noexcept
    : lower(grid_lower), upper(grid_upper), guard(grid_guard),
      delta( (grid_upper - grid_lower) / grid_N ) {}

  // abscissa
  constexpr T absc( int i, T shift_from_lb ) noexcept {
    return lower + delta * ( i + shift_from_lb - guard );
  }

  // constexpr T rel_absc( int i, T shift_from_lb ) noexcept {
  //   return (lower / delta) - guard + i + shift_from_lb;
  // }

  // constexpr int index ( T abscissa ) noexcept {
  //   return static_cast<int>( ( abscissa - lower ) / delta + guard );
  // }

  // constexpr int index_from_rel ( T rel_abscissa ) noexcept {
  //   return static_cast<int>( rel_abscissa - (lower / delta) + guard );
  // }
};

// Grid is designed to represent the supergrid
template < std::size_t DGrid, typename T = Real  >
using Grid = std::array< Gridline<T>, DGrid >;

#include "vector.hpp"
namespace mem {
  template < class Grid_t  >
  vec_def_member_getter( const Grid_t, delta );

  template < class Grid_t  >
  vec_def_member_getter( const Grid_t, lower );

  template < class Grid_t  >
  vec_def_member_getter( const Grid_t, upper );

  template < class Grid_t  >
  vec_def_member_getter( const Grid_t, guard );
}


#endif
