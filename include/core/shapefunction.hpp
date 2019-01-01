#ifndef _SHAPEFUNCTION_HPP_
#define _SHAPEFUNCTION_HPP_

struct ShapeFunction {
  constexpr int support() const {
    if constexpr ( traits::shape == "Nearest Grid Point" ) return 1;
    else if ( traits::shape == "Cloud In Cell" ) return 2;
    else if ( traits::shape == "Triangular Cloud" ) return 3;
    else if ( traits::shape == "Piecewise Cubic Spline" ) return 4;
    // TODOL add check for unknown type. Note static_assert cannot be used here
    // else static_assert(false, "Unknown Particle Shape");
  }

  constexpr Real radius() const {
    return support() * 0.5;
  }

  constexpr Real weight( Real dx ) const {
    // TODOL check correctness Real&&
    Real&& absdx = std::abs( std::move(dx) );

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
