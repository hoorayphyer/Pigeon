#ifndef  _KNL_COORDINATE_HPP_
#define  _KNL_COORDINATE_HPP_

#include "kernel/coordsys_predef.hpp"
#include "apt/numeric.hpp"
#include "apt/vec.hpp"
#include "apt/virtual_vec.hpp"

inline constexpr
long double PI_CONST =3.141592653589793238462643383279502884197169399375105820974944592307816406286L;

namespace knl {

  template < coordsys Coord = coordsys::Cartesian >
  struct coord {
    template < int N, typename T >
    static inline T h ( T x1, T x2, T x3 ) {
      static_assert( N == 1 || N == 2 || N == 3, "h index out of bounds" );
      return 1.0;
    }

    template < int NN, typename T >
    static inline T hh ( T x1, T x2, T x3 ) {
      static_assert( NN == 12 || NN == 23 || NN == 31, "hh index out of bounds");
      return 1.0;
    }

    template < typename T >
    static inline T hhh ( T x1, T x2, T x3 ) { return 1.0; }

    template < class X, class V, typename T >
    static inline auto geodesic_move( X& x, const apt::VecExpression<V>& v, const T& dt ) {
      apt::Vec<decltype(v[0] * dt), apt::ndim_v<V>> dx = v * dt;
      x += dx;
      return dx;
    }
  };

  template <>
  struct coord < coordsys::LogSpherical > {
    template < int N, typename T >
    static inline T h ( T logr, T theta, T phi ) {
      static_assert( N == 1 || N == 2 || N == 3, "h index out of bounds" );

      if constexpr ( N == 3 ) return std::exp(logr) * std::sin(theta);
      else return std::exp(logr);
    }

    template < int NN, typename T >
    static inline T hh ( T logr, T theta, T phi ) {
      static_assert( NN == 12 || NN == 23 || NN == 31, "hh index out of bounds");

      if constexpr ( NN == 12 ) return std::exp( 2.0 * logr);
      else return std::exp( 2.0 * logr) * std::sin(theta);
    }

    template < typename T >
    static inline T hhh ( T logr, T theta, T phi ) {
      return std::exp( 3.0 * logr ) * std::sin(theta);
    }

    template < class X, class V, typename T >
    static inline auto geodesic_move( X& x, V& v, const T& dt ) {
      // TODOL: in implementing this, we assumed 1) no crossing through center and 2) no crossing through symmetry axes

      constexpr T PI = PI_CONST;
      // dx is the return value and meanwhile it serves as a temporary
      apt::Vec<T, apt::ndim_v<V>> dx = v * dt;

      auto& dlogr = dx[0];
      auto& dtheta = dx[1];
      auto& dphi = dx[2];

      const auto& logr  = x[0];
      const auto& theta = x[1];

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

      { // compute dx in this block
        rot.set_angle(theta);

        dlogr += std::exp( logr ); // now dx is ( r0 + v_r * dt, v_theta * dt, v_phi * dt )

        auto r_final = apt::abs( dx );

        dtheta = std::acos( ( dlogr * rot.cos() - dtheta * rot.sin() ) / r_final ) - theta;
        dphi = std::atan( dphi / std::abs( dlogr * rot.sin() + dtheta * rot.cos() ) );

        dlogr = std::log( r_final ) - logr;
      }
      { // rebase velocity to the new position
        // first rotate in r-theta plane to equator. Location is rotated by PI/2 - theta, so velocity components are rotated by theta - PI/2.
        rot.set_angle( theta - PI / 2.0 );
        rot.rotate( v[0], v[1] );

        // then rotate in r-phi plane. Location is rotated by dphi, so velocity components are rotated -dphi
        rot.set_angle( -dphi );
        rot.rotate( v[0], v[2] );

        // last rotate in r-theta plane to new location. Location is rotated by theta_new - PI/2, so velocity components are rotated by PI/2 - theta_new
        rot.set_angle( PI / 2.0 - theta - dtheta );
        rot.rotate( v[0], v[1] );
      }

      { // in this block we update x
        // NOTE we don't normalize cyclic coordinates such as \phi in spherical. This operation is designated to mesh_shape_interplay
        x += dx;
      }

      return dx;
    }
  };
}




#endif
