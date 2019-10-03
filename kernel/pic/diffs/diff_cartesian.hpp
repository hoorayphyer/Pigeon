#ifndef _FIELD_DERIVATIVE_CARTESIAN_HPP_
#define _FIELD_DERIVATIVE_CARTESIAN_HPP_

#include "field/field.hpp"
#include "manifold/grid.hpp"
#include "field/yee.hpp"

namespace field {
  template < int DGrid, typename T >
  constexpr T diff_one( T , T , T ) noexcept {
    return 1.0;
  }

  // NOTE all derivatives are second order
  template < int DGrid, typename T, int I, offset_t f_ofs >
  constexpr T diff( const T& f, T , T , const mani::Grid<T,DGrid>& g, const apt::Index<DGrid>& s ) noexcept {
    static_assert( DGrid > 1 );
    static_assert( I >= 0 && I < 3 );
    if constexpr ( I >= DGrid ) return 0.0;

    if constexpr ( INSITU == f_ofs ) return ( *(&f + s[I]) - f ) / g[I].delta();
    else return ( f - *(&f - s[I]) ) / g[I].delta();
  }

  template < int DGrid, typename T, offset_t Ftype >
  constexpr T (*Diff( int Fcomp, int Icoord )) (const T& f, T, T, const mani::Grid<T,DGrid>& g, const apt::Index<DGrid>& s)  noexcept {
    if ( 0 == Fcomp ) {
      if ( 1 == Icoord ) return diff<DGrid,T,1,!Ftype>;
      if ( 2 == Icoord ) return diff<DGrid,T,2,!Ftype>;
    } else if ( 1 == Fcomp ) {
      if ( 2 == Icoord ) return diff<DGrid,T,2,!Ftype>;
      if ( 0 == Icoord ) return diff<DGrid,T,0,!Ftype>;
    } else {
      if ( 0 == Icoord ) return diff<DGrid,T,0,!Ftype>;
      if ( 1 == Icoord ) return diff<DGrid,T,1,!Ftype>;
    }
    return nullptr;
  }
}

#endif
