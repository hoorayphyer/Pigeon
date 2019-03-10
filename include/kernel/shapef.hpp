#ifndef _KNL_SHAPEF_HPP_
#define _KNL_SHAPEF_HPP_

#include "kernel/shape_predef.hpp"

namespace knl {
  template < shape S >
  struct shapef_t {
    static const int support = static_cast<int>(S);

    template < typename T >
    constexpr T operator() ( T dx ) const noexcept {
      dx = std::abs(dx);

      if constexpr( shape::Nearest_Grid_Point == S ) {
          return static_cast<T> ( dx < 0.5 );
        }
      else if ( shape::Cloud_In_Cell == S ) {
        return std::max ( 1.0 - dx, 0.0 );
      }
      else if ( shape::Triangular_Cloud == S ) {
        return
          static_cast<T>( dx < 0.5 ) * ( 0.75 - dx * dx )
          +
          static_cast<T>( dx >= 0.5 && dx < 1.5 ) * 0.5 * ( 1.5 - dx ) * ( 1.5 - dx );
      }
      else if ( shape::Piecewise_Cubic_Spline == S ) {
        return
          static_cast<T>( dx < 1.0 ) * ( 2.0 / 3.0 - dx * dx * ( 1.0 - 0.5 * dx ) )
          +
          static_cast<T>( dx >= 1.0 && dx < 2.0 ) * ( 2.0 - dx ) * ( 2.0 - dx ) * ( 2.0 - dx )  / 6.0;
      }
    }

  };

}


#endif
