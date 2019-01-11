#ifndef _SHAPEFUNCTION_HPP_
#define _SHAPEFUNCTION_HPP_

#include "types.hpp"

template < ct_string shape, typename T = Real >
struct ShapeFunction {
  static constexpr int support() {
    if constexpr ( shape == "Nearest Grid Point" ) return 1;
    else if ( shape == "Cloud In Cell" ) return 2;
    else if ( shape == "Triangular Cloud" ) return 3;
    else if ( shape == "Piecewise Cubic Spline" ) return 4;
    // TODOL add check for unknown type. Note static_assert cannot be used here
    else static_assert(false, "Unknown Particle Shape");
  }

  constexpr T radius = support() * 0.5;

  constexpr T operator() ( T dx ) {
    dx = std::abs(dx);

    if constexpr( shape == "Nearest_Grid_Point" ) {
        return static_cast<T> ( dx <= 0.5 );
      }

    else if ( shape == "Cloud In Cell" ) {
        return std::max ( 1.0 - dx, 0.0 );
      }


    else if ( shape == "Triangular_Cloud" ) {
        return
          static_cast<T>( dx < 0.5 ) * ( 0.75 - dx * dx )
          +
          static_cast<T>( dx >= 0.5 && dx < 1.5 ) * 0.5 * ( 1.5 - dx ) * ( 1.5 - dx );
      }

    else if ( shape == "Piecewise_Cubic_Spline" ) {
        return
          static_cast<T>( dx < 1.0 ) * ( 2.0 / 3.0 - dx * dx * ( 1.0 - 0.5 * dx ) )
          +
          static_cast<T>( dx >= 1.0 && dx < 2.0 ) * ( 2.0 - dx ) * ( 2.0 - dx ) * ( 2.0 - dx )  / 6.0;
      }

  }

};

#endif
