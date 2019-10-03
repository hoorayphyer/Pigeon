#ifndef  _MANIFOLD_CURVILINEAR_IMPL_HPP_
#define  _MANIFOLD_CURVILINEAR_IMPL_HPP_

#include "manifold/curvilinear.hpp"
#include "apt/numeric.hpp"
#include "apt/vec.hpp"
#include "apt/virtual_vec.hpp"

template < typename T >
inline constexpr T PI = std::acos(-1.0L);

namespace mani {
  struct CartesianCoordSys {
    template < int N, typename T >
    static inline T h ( T = 0.0, T = 0.0, T = 0.0 ) noexcept {
      return 1.0;
    }

    template < int N, typename T >
    static inline T hh ( T = 0.0, T = 0.0, T = 0.0 ) noexcept {
      static_assert( N == 0 || N == 1 || N == 2, "h index out of bounds" );
      return 1.0;
    }

    template < typename T >
    static inline T hhh ( T = 0.0, T = 0.0, T = 0.0 ) noexcept { return 1.0; }

    template < class X, class P, typename T >
    static inline auto geodesic_move( X& x, P& p, T dt, bool is_massive ) noexcept {
      apt::array<T, P::NDim> dx;
      dt /= std::sqrt( is_massive + apt::sqabs(p) );
      for ( int i = 0; i < P::NDim; ++i ) {
        dx[i] = p[i] * dt;
        x[i] += dx[i];
      }
      return dx;
    }
  };

  struct LogSphericalCoordSys {
    template < int N, typename T >
    static inline T h ( T logr, T theta, T = 0.0 ) noexcept {
      static_assert( N == 0 || N == 1 || N == 2, "h index out of bounds" );
      if constexpr ( N == 2 ) return std::exp(logr) * std::sin(theta);
      else return std::exp(logr);
    }

    template < int N, typename T >
    static inline T hh ( T logr, T theta, T = 0.0 ) noexcept {
      static_assert( N == 0 || N == 1 || N == 2, "h index out of bounds" );
      if constexpr ( N == 2 ) return std::exp( 2.0 * logr);
      else return std::exp( 2.0 * logr) * std::sin(theta);
    }

    template < typename T >
    static inline T hhh ( T logr, T theta, T = 0.0 ) noexcept {
      return std::exp( 3.0 * logr ) * std::sin(theta);
    }

    template < class X, class P, typename T >
    static inline auto geodesic_move( X& x, P& p, T dt, bool is_massive ) noexcept {
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

        dt = (dt == static_cast<T>(0)) ? 0.0 : std::atan( dx[2] / dt ) + PI<T> * is_massive;
        dx[2] = std::sqrt( dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
        dx[1] = std::acos( dx[0] / dx[2] );
        std::swap( dt, dx[2] );
        dt = std::log(dt);
        // by now, dt = ln(r_final), dx[0] is free, dx[1] = \theta_final, dx[2] = \phi_final - \phi_init
      }
      { // rebase velocity to the new position
        // first rotate in r-theta plane to equator. Location is rotated by PI/2 - theta, so velocity components are rotated by theta - PI/2.
        T pmg = apt::abs(p); // explicitly conserve energy
        cos = std::cos( x[1] - PI<T> / 2.0 );
        sin = std::sin( x[1] - PI<T> / 2.0 );
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
        cos = std::cos( PI<T> / 2.0 - dx[1] );
        sin = std::sin( PI<T> / 2.0 - dx[1] );
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
        dx[1] += is_massive ? ( -x[1] + ( x[1] > PI<T> / 2.0 ) * 2.0 * PI<T> ) : x[1] ;

        x[2] += dx[2];
        dx[2] -= PI<T> * is_massive;

      }

      return dx;
    }
  };
}




#endif
