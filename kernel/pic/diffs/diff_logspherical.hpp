#ifndef _FIELD_DERIVATIVE_LOGSPHERICAL_HPP_
#define _FIELD_DERIVATIVE_LOGSPHERICAL_HPP_

#include "manifold/grid.hpp"
#include "manifold/curvilinear_impl.hpp"

using Metric = mani::LogSphericalCoordSys;

namespace field {
  template < int DGrid, typename T >
  constexpr T diff_zero( const T& f, T lnr, T theta, const mani::Grid<T,DGrid>& g, const apt::Index<DGrid+1>& s ) noexcept {
    return 0.0;
  }
}

template < int DGrid, typename T >
constexpr T diff_one( T , T , T ) noexcept {
  return 1.0;
}

namespace field {
  template < int DGrid, typename T, int Fcomp, int I >
  struct FH {
    static_assert( DGrid > 1 && I >= 0 && I < 3 );

    const T& f;
    const T& lnr;
    const T& theta;
    const mani::Grid<T,DGrid>& g;
    const apt::Index<DGrid+1>& s;

    constexpr T operator()(int shift) const noexcept {
      return *(&f + shift * s[I]) * Metric::template h<Fcomp>(lnr + (0==I) * shift * g[I].delta(), theta + (1==I) * shift * g[I].delta() );
    }
  };

  // NOTE all derivatives are second order
  template < int DGrid, typename T, int Fcomp, int I, offset_t f_ofs >
  constexpr T diff( const T& f, T lnr, T theta, const mani::Grid<T,DGrid>& g, const apt::Index<DGrid+1>& s ) noexcept {
    if constexpr ( I >= DGrid ) return 0.0;

    FH<DGrid,T,Fcomp,I> fh{f,lnr,theta,g,s};
    if constexpr ( yee::Btype == f_ofs )
      return (fh(1) - fh(0)) / g[I].delta();
    else
      return (fh(0) - fh(-1)) / g[I].delta();
  }

  template < int DGrid, typename T, offset_t Ftype >
  constexpr T (*Diff( int Fcomp, int Icoord )) (const T& f, T lnr, T theta, const mani::Grid<T,DGrid>& g, const apt::Index<DGrid+1>& s)  noexcept {
    if ( 0 == Fcomp ) {
      if ( 1 == Icoord ) return diff<DGrid,T,0,1,Ftype>;
      if ( 2 == Icoord ) return diff<DGrid,T,0,2,Ftype>;
    } else if ( 1 == Fcomp ) {
      if ( 2 == Icoord ) return diff<DGrid,T,1,2,Ftype>;
      if ( 0 == Icoord ) return diff<DGrid,T,1,0,Ftype>;
    } else {
      if ( 0 == Icoord ) return diff<DGrid,T,2,0,Ftype>;
      if ( 1 == Icoord ) return diff<DGrid,T,2,1,Ftype>;
    }
    return nullptr;
  }

  template < typename T >
  constexpr T(*HHk(int i) )( T,T,T ) noexcept {
    switch (i) {
    case 0 : return Metric::hh<0,T>;
    case 1 : return Metric::hh<1,T>;
    case 2 : return Metric::hh<2,T>;
    }
    return nullptr;
  }

  template < typename T >
  constexpr T(*HiHj(int i, int j) )( T,T,T ) noexcept {
    return HHk<T>(3-i-j);
  }

  template < int DGrid, typename T, int Fcomp, int I >
  constexpr T diff_1sided_from_right( const T& f, T lnr, T theta, const mani::Grid<T,DGrid>& g, const apt::Index<DGrid+1>& s ) noexcept {
    if constexpr ( I >= DGrid ) return 0.0;
    FH<DGrid,T,Fcomp,I> fh{f,lnr,theta,g,s};
    return ( - 2.0 * fh(0) + 3.0 * fh(1) - fh(2) ) / g[I].delta();
  }

  template < int DGrid, typename T, int Fcomp, int I >
  constexpr T diff_1sided_from_left( const T& f, T lnr, T theta, const mani::Grid<T,DGrid>& g, const apt::Index<DGrid+1>& s ) noexcept {
    if constexpr ( I >= DGrid ) return 0.0;
    // NOTE For B, the upper boundary is actually a guard cell. We use the same convention that f is in the same cell as the evaluated derivative
    FH<DGrid,T,Fcomp,I> fh{f,lnr,theta,g,s};
    return ( 2.0 * fh(-1) - 3.0 * fh(-2) + fh(2) ) / g[I].delta();
  }

  template < int DGrid, typename T, bool IsUpper >
  constexpr T diff_axis_Ephi_theta( const T& f, T lnr, T theta, const mani::Grid<T,DGrid>& g, const apt::Index<DGrid+1>& s ) noexcept {
    constexpr int Fcomp = 2;
    constexpr int I = 1;
    if constexpr ( I >= DGrid ) return 0.0;

    FH<DGrid,T,Fcomp,I> fh{f,lnr,theta,g,s};
    if constexpr (IsUpper) return ( - 81 * fh(-1) + fh(-2) ) / std::pow( 3*std::exp(lnr)*g[I].delta(), 2.0);
    else return ( 81 * fh(0) - fh(1) ) / std::pow( 3*std::exp(lnr)*g[I].delta(), 2.0);
  }

  template < int DGrid, typename T, int Fcomp, int I, offset_t f_ofs >
  constexpr T diff_beyond_hi_axis( const T& f, T lnr, T theta, const mani::Grid<T,DGrid>& g, const apt::Index<DGrid+1>& s ) noexcept {
    constexpr bool IsAxissymmetric = ( Fcomp + I == 3 );
    return ( IsAxissymmetric ? 1 : -1 ) * diff<DGrid,T,Fcomp,I,f_ofs>( *(&f-s[1]), lnr, theta - g[1].delta(), g, s );
  }

  template < int DGrid, typename T, int Fcomp, int I >
  constexpr T diff_1sided_from_right_beyond_hi_axis( const T& f, T lnr, T theta, const mani::Grid<T,DGrid>& g, const apt::Index<DGrid+1>& s ) noexcept {
    constexpr bool IsAxissymmetric = ( Fcomp + I == 3 );
    return ( IsAxissymmetric ? 1 : -1 ) * diff_1sided_from_right<DGrid,T,Fcomp,I>( *(&f-s[1]), lnr, theta - g[1].delta(), g, s );
  }

  template < int I, typename T >
  constexpr T hh_beyond_hi_axis( T lnr, T theta, T = 0 ) noexcept {
    return std::abs(Metric::template hh<I,T>(lnr, theta, 0));
  }

}

#endif
