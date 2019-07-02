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
    static inline T h ( T = 0.0, T = 0.0, T = 0.0 ) {
      return 1.0;
    }

    template < int N, typename T >
    static inline T hh ( T = 0.0, T = 0.0, T = 0.0 ) {
      static_assert( N == 0 || N == 1 || N == 2, "h index out of bounds" );
      return 1.0;
    }

    template < typename T >
    static inline T hhh ( T = 0.0, T = 0.0, T = 0.0 ) { return 1.0; }

    template < class X, class P, typename T >
    static inline auto geodesic_move( X& x, const apt::VecExpression<P>& p, T dt, bool is_massive ) noexcept {
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
    static inline T h ( T logr, T theta, T = 0.0 ) {
      static_assert( N == 0 || N == 1 || N == 2, "h index out of bounds" );
      if constexpr ( N == 2 ) return std::exp(logr) * std::sin(theta);
      else return std::exp(logr);
    }

    template < int N, typename T >
    static inline T hh ( T logr, T theta, T = 0.0 ) {
      static_assert( N == 0 || N == 1 || N == 2, "h index out of bounds" );
      if constexpr ( N == 2 ) return std::exp( 2.0 * logr);
      else return std::exp( 2.0 * logr) * std::sin(theta);
    }

    template < typename T >
    static inline T hhh ( T logr, T theta, T = 0.0 ) {
      return std::exp( 3.0 * logr ) * std::sin(theta);
    }

    template < class X, class P, typename T >
    static inline auto geodesic_move( X& x, P& p, T dt, bool is_massive ) noexcept {
      apt::array<T, P::NDim> dx;

      // TODOL: in implementing this, we assumed 1) no crossing through center and 2) no crossing through symmetry axes
      class Rotator {
      private:
        T _cos = 0.0;
        T _sin = 0.0;
      public:
        inline void set_angle( T angle ) noexcept {
          _cos = std::cos(angle);
          _sin = std::sin(angle);
        }

        // in-place rotation
        inline void rotate( T& v1, T& v2 ) const noexcept {
          v1 = _cos * v1 - _sin * v2;
          v2 = ( _sin * v1 + v2 ) / _cos;
        }

        inline T cos() const noexcept { return _cos; }
        inline T sin() const noexcept { return _sin; }
      } rot;

      dt /= std::sqrt( is_massive + apt::sqabs(p) );
      for ( int i = 0; i < P::NDim; ++i )
        dx[i] = p[i] * dt;
      // now dx holds displacements under the local cartesian frame

      { // compute final coordiantes. In this section, x remains untouched
        rot.set_angle(x[1]);
        dx[0] += std::exp( x[0] ); // dx[0] = r_i + d_r
        dt = dx[0] * rot.sin() + dx[1] * rot.cos(); // NOTE now dt stores an important determinant
        is_massive = ( dt < 0 ); // NOTE is_massive now stores whether axis crossing happened when displacing in the theta direction
        dt = std::atan( dx[2] / dt ) + PI<T> * is_massive;
        dx[2] = std::sqrt( dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);

        dx[0] *= rot.cos();
        dx[1] *= rot.sin();

        dx[1] = ( dx[0]*dx[0] - dx[1]*dx[1] ) / ( dx[0] + dx[1] );
        dx[1] = std::acos( dx[1] / dx[2] );
        std::swap( dt, dx[2] );
        dt = std::log(dt);
        // by now, dt = ln(r_final), dx[0] is free, dx[1] = \theta_final, dx[2] = \phi_final - \phi_init
      }
      { // rebase velocity to the new position
        // first rotate in r-theta plane to equator. Location is rotated by PI/2 - theta, so velocity components are rotated by theta - PI/2.
        rot.set_angle( x[1] - PI<T> / 2.0 );
        rot.rotate( p[0], p[1] );

        // then rotate in r-phi plane. Location is rotated by dphi, so velocity components are rotated -dphi
        rot.set_angle( -dx[2] );
        rot.rotate( p[0], p[2] );

        // last rotate in r-theta plane to new location. Location is rotated by theta_new - PI/2, so velocity components are rotated by PI/2 - theta_new
        rot.set_angle( PI<T> / 2.0 - dx[1] );
        rot.rotate( p[0], p[1] );
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
