#ifndef  _GRID_HPP_
#define  _GRID_HPP_

// grid here is the supergrid in one dimension
struct Grid {
  Real lower;
  Real upper;
  int guard;
  Real delta;

  constexpr Grid() : Grid( 0.0, 1.0, 0, 1 ) {}

  constexpr Grid( Real grid_lower, Real grid_upper, int grid_guard, int grid_N )
    : lower(grid_lower), upper(grid_upper), guard(grid_guard),
      delta( (grid_upper - grid_lower) / grid_N ) {}

  constexpr Real abscissa( int i, Real shift_from_lb ) noexcept {
    return ( lower - delta * guard + shift_from_lb * delta ) + i * delta;
  }
};


#endif
