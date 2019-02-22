#ifndef  _KNL_COORDINATE_HPP_
#define  _KNL_COORDINATE_HPP_

#include "kernel/coordsys_predef.hpp"
#include "apt/numeric.hpp"
#include "apt/vec.hpp"

template < typename T >
inline constexpr T PI (3.141592653589793238462643383279502884197169399375105820974944592307816406286L);

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
    static inline auto geodesic_move( apt::VecExpression<X>& x, const apt::VecExpression<V>& v, const T& dt ) {
      apt::Vec<decltype(std::get<0>(v) * dt), V::size> dx = v * dt;
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
    static inline auto geodesic_move( apt::VecExpression<X>& x, apt::VecExpression<V>& v, const T& dt ) {
      // TODOL: in implementing this, we assumed 1) no crossing through center and 2) no crossing through symmetry axes

      // dx is the return value and meanwhile it serves as a temporary
       apt::Vec<decltype(std::get<0>(v) * dt), V::size> dx = v * dt;
      auto& dlogr = std::get<0>(dx);
      auto& dtheta = std::get<1>(dx);
      auto& dphi = std::get<2>(dx);

      const auto& logr  = std::get<0>(x);
      const auto& theta = std::get<1>(x);
      { // compute dx in this block
        dlogr += std::exp( logr ); // after this step, dx should be ( r0 + v_r * dt, v_theta * dt, v_phi * dt )

        auto r_final = apt::abs( dx );

        dtheta = std::acos( ( dlogr * std::cos(theta) - dtheta * std::sin(theta) ) / r_final ) - theta;
        dphi = std::atan( dphi / std::abs( dlogr * std::sin(theta) + dtheta * std::cos(theta) ) );

        dlogr = std::log(r_final) - logr;
      }
      { // rebase velocity to the new position
        auto cdp = std::cos(dphi);
        auto sdp = std::sin(dphi);
        auto sT = std::sin( 2*theta + dtheta );
        auto cT = std::cos( 2*theta + dtheta );
        auto sdt = std::sin( dtheta );
        auto cdt = std::cos( dtheta );

        // store old values of v
        auto v_r = std::get<0>(v);
        auto v_theta = std::get<1>(v);
        auto v_phi = std::get<2>(v);

        std::get<0>(v) =
          v_r * 0.5 * ( ( 1 - cdp ) * cT + ( 1 + cdp )  * cdt )
          +  v_theta  * 0.5 * ( ( cdp - 1 ) * sT+ ( 1 + cdp ) * sdt )
          + v_phi * std::sin( theta + dtheta ) * sdp;

        std::get<1>(v) =
          v_r * 0.5 * ( ( cdp - 1 ) * sT - ( 1 + cdp ) * sdt )
          + v_theta * 0.5 * ( ( cdp - 1 ) * cT + ( 1 + cdp ) * cdt )
          + v_phi * std::cos( theta + dtheta ) * sdp;

        std::get<2>(v) =
          - v_r * std::sin(theta) * sdp
          - v_theta * std::cos(theta) * sdp
          + v_phi * cdp;
      }

      { // in this block we update x
        x += dx;
        // normalize phi
        auto& phi = std::get<2>(x);
        phi -= 2 * PI<T> * std::floor( phi / (2 * PI<T>) );
      }

      return dx;
    }
  };
}




#endif
