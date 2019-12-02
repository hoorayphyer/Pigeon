#ifndef _METRIC_LOGSPHERICAL_HPP_
#define _METRIC_LOGSPHERICAL_HPP_

#include <cmath>

namespace metric {
  template < typename T >
  struct LogSpherical {
    static constexpr T PI = std::acos(-1.0L);

    template < int N >
    static inline T h ( T logr, T theta, T = 0.0 ) noexcept {
      static_assert( N == 0 || N == 1 || N == 2, "h index out of bounds" );
      if constexpr ( N == 2 ) return std::exp(logr) * std::sin(theta);
      else return std::exp(logr);
    }

    template < int N >
    static inline T hh ( T logr, T theta, T = 0.0 ) noexcept {
      static_assert( N == 0 || N == 1 || N == 2, "h index out of bounds" );
      if constexpr ( N == 2 ) return std::exp( 2.0 * logr);
      else return std::exp( 2.0 * logr) * std::sin(theta);
    }

    static inline T hhh ( T logr, T theta, T = 0.0 ) noexcept {
      return std::exp( 3.0 * logr ) * std::sin(theta);
    }

    template < class X, class P >
    static auto geodesic_move( X& x, P& p, T dt, bool is_massive ) noexcept;
  };
}

#include "apt/index.hpp"
#include "apt/grid.hpp"

namespace field {
  template < int DGrid, typename T >
  constexpr T diff_zero( const T& f, T lnr, T theta, const apt::Grid<T,DGrid>& g, const apt::Index<DGrid+1>& s ) noexcept {
    return 0.0;
  }

  template < typename T >
  constexpr T diff_one( T , T , T ) noexcept {
    return 1.0;
  }
}


#include "field/yee.hpp"
#include "field/field.hpp"

namespace field {
  template < int DGrid, typename T, int Fcomp, int I >
  struct FH {
    using Metric = metric::LogSpherical<T>;
    static_assert( DGrid > 1 && I >= 0 && I < 3 );

    const T& f;
    const T& lnr;
    const T& theta;
    const apt::Grid<T,DGrid>& g;
    const apt::Index<DGrid+1>& s;

    constexpr T operator()(int shift) const noexcept {
      return *(&f + shift * s[I]) * Metric::template h<Fcomp>(lnr + (0==I) * shift * g[0].delta(), theta + (1==I) * shift * g[1].delta() );
    }
  };

  // NOTE all derivatives are second order
  template < int DGrid, typename T, int Fcomp, int I, offset_t f_ofs >
  constexpr T diff( const T& f, T lnr, T theta, const apt::Grid<T,DGrid>& g, const apt::Index<DGrid+1>& s ) noexcept {
    if constexpr ( I >= DGrid ) return 0.0;

    FH<DGrid,T,Fcomp,I> fh{f,lnr,theta,g,s};
    if constexpr ( yee::Btype == f_ofs )
      return (fh(1) - fh(0)) / g[I].delta();
    else
      return (fh(0) - fh(-1)) / g[I].delta();
  }

  template < int DGrid, typename T, offset_t Ftype >
  constexpr T (*Diff( int Fcomp, int Icoord )) (const T& f, T lnr, T theta, const apt::Grid<T,DGrid>& g, const apt::Index<DGrid+1>& s)  noexcept {
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
    using Metric = metric::LogSpherical<T>;
    switch (i) {
    case 0 : return Metric::template hh<0>;
    case 1 : return Metric::template hh<1>;
    case 2 : return Metric::template hh<2>;
    }
    return nullptr;
  }

  template < typename T >
  constexpr T(*HiHj(int i, int j) )( T,T,T ) noexcept {
    return HHk<T>(3-i-j);
  }

  template < int DGrid, typename T, int Fcomp, int I >
  constexpr T diff_1sided_from_right( const T& f, T lnr, T theta, const apt::Grid<T,DGrid>& g, const apt::Index<DGrid+1>& s ) noexcept {
    if constexpr ( I >= DGrid ) return 0.0;
    FH<DGrid,T,Fcomp,I> fh{f,lnr,theta,g,s};
    return ( - 2.0 * fh(0) + 3.0 * fh(1) - fh(2) ) / g[I].delta();
  }

  template < int DGrid, typename T, bool IsUpper >
  constexpr T diff_axis_Ephi_theta( const T& f, T lnr, T theta, const apt::Grid<T,DGrid>& g, const apt::Index<DGrid+1>& s ) noexcept {
    constexpr int Fcomp = 2;
    constexpr int I = 1;
    if constexpr ( I >= DGrid ) return 0.0;

    if constexpr (IsUpper) return 2 * (*( &f - s[I] )) * std::exp(-lnr) / std::sin(theta); // NOTE theta here is larger than pi, and this formula is correct.
    else return 2 * f * std::exp(-lnr) / std::sin(theta);
  }

}


#include "apt/numeric.hpp"
#include "apt/vec.hpp"
#include "apt/virtual_vec.hpp"

namespace metric {
  template < typename T >
  template < class X, class P >
  auto LogSpherical<T>::geodesic_move( X& x, P& p, T dt, bool is_massive ) noexcept {
    apt::array<T, P::NDim> dx;

    // TODOL: in implementing this, we assumed no crossing through center

    dt /= std::sqrt( is_massive + apt::sqabs(p) );
    for ( int i = 0; i < P::NDim; ++i )
      dx[i] = p[i] * dt;
    // now dx holds displacements under the local cartesian frame

    T cos = 0.0;
    T sin = 0.0;
    { // compute final coordiantes. In this section, x remains untouched
      cos = std::cos(x[1]);
      sin = std::sin(x[1]);
      dx[0] += std::exp( x[0] ); // dx[0] = r_i + d_r
      dt = dx[0] * sin + dx[1] * cos; // NOTE now dt stores an important determinant
      is_massive = ( dt < 0 ); // NOTE is_massive now stores whether axis crossing happened when displacing in the theta direction

      // rotate dx[0], dx[1] by x[1]. This solves the loss of significance under some circumstances when compiler optimization is strong, which results in either std::sqrt-ing a negative number or std::acos-ing a number larger than 1, causing NAN error
      dx[0] = dx[0] * cos - dx[1] * sin;
      dx[1] = dt;

      dt = (dt == static_cast<T>(0)) ? 0.0 : std::atan( dx[2] / dt ) + PI * is_massive;
      dx[2] = std::sqrt( dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
      dx[1] = std::acos( dx[0] / dx[2] );
      std::swap( dt, dx[2] );
      dt = std::log(dt);
      // by now, dt = ln(r_final), dx[0] is free, dx[1] = \theta_final, dx[2] = \phi_final - \phi_init
    }
    { // rebase velocity to the new position
      // first rotate in r-theta plane to equator. Location is rotated by PI/2 - theta, so velocity components are rotated by theta - PI/2.
      T pmg = apt::abs(p); // explicitly conserve energy
      cos = std::cos( x[1] - PI / 2.0 );
      sin = std::sin( x[1] - PI / 2.0 );
      dx[0] = p[0];
      p[0] = dx[0] * cos - p[1] * sin;
      p[1] = dx[0] * sin + p[1] * cos;

      // then rotate in r-phi plane. Location is rotated by dphi, so velocity components are rotated -dphi
      cos = std::cos( -dx[2] );
      sin = std::sin( -dx[2] );
      dx[0] = p[0];
      p[0] = dx[0] * cos - p[2] * sin;
      p[2] = dx[0] * sin + p[2] * cos;

      // last rotate in r-theta plane to new location. Location is rotated by theta_new - PI/2, so velocity components are rotated by PI/2 - theta_new
      cos = std::cos( PI / 2.0 - dx[1] );
      sin = std::sin( PI / 2.0 - dx[1] );
      dx[0] = p[0];
      p[0] = dx[0] * cos - p[1] * sin;
      p[1] = dx[0] * sin + p[1] * cos;
      p *= ( pmg / apt::abs(p) );
    }

    { // in this block we update x
      // NOTE we don't normalize cyclic coordinates such as \phi in spherical. This operation is designated to mesh_shape_interplay
      dx[0] = dt - x[0];
      std::swap( dt, x[0] );

      std::swap( x[1], dx[1] );
      dx[1] *= -1;
      dx[1] += is_massive ? ( -x[1] + ( x[1] > PI / 2.0 ) * 2.0 * PI ) : x[1] ;

      x[2] += dx[2];
      dx[2] -= PI * is_massive;

    }

    return dx;
  }
}


#endif
