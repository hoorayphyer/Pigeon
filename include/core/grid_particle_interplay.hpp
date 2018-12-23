#ifndef _GRIDPARTICLEINTERPLAY_HPP_
#define _GRIDPARTICLEINTERPLAY_HPP_

#include "grid.hpp"

template < int I_Dim >
struct shape_func {
private:
public:
  static_assert( I_Dim < traits::Dgrid, "dimension overflow");

  static constexpr int support() {
    if constexpr ( traits::shape == "Nearest Grid Point" ) return 1;
    else if ( traits::shape == "Cloud In Cell" ) return 2;
    else if ( traits::shape == "Triangular Cloud" ) return 3;
    else if ( traits::shape == "Piecewise Cubic Spline" ) return 4;
    // TODO add check for unknown type. Note static_assert cannot be used here
    // else static_assert(false, "Unknown Particle Shape");
  }

  static constexpr Real radius() {
    return static_cast<Real>( support() ) * grid<I_Dim>::delta * 0.5;
  }

  static constexpr Real weight( const Real delta_q ) {
    const Real&& absdx = std::abs( std::move(delta_q) / grid<I_Dim>::delta );

    if constexpr( traits::shape == "Nearest_Grid_Point" ) {
        return static_cast<Real> ( absdx <= 0.5 );
      }

    else if ( traits::shape == "Cloud In Cell" ) {
        return std::max ( 1.0 - absdx, 0.0 );
      }


    else if ( traits::shape == "Triangular_Cloud" ) {
        return
          static_cast<Real>( absdx < 0.5 ) * ( 0.75 - absdx * absdx )
          +
          static_cast<Real>( absdx >= 0.5 && absdx < 1.5 ) * 0.5 * ( 1.5 - absdx ) * ( 1.5 - absdx );
      }

    else if ( traits::shape == "Piecewise_Cubic_Spline" ) {
        return
          static_cast<Real>( absdx < 1.0 ) * ( 2.0 / 3.0 - absdx * absdx * ( 1.0 - 0.5 * absdx ) )
          +
          static_cast<Real>( absdx >= 1.0 && absdx < 2.0 ) * ( 2.0 - absdx ) * ( 2.0 - absdx ) * ( 2.0 - absdx )  / 6.0;
      }

  }

};



#endif
