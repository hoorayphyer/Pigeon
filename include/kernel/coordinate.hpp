#ifndef  _KNL_COORDINATE_HPP_
#define  _KNL_COORDINATE_HPP_

#include "apt/numeric.hpp"

template < typename T >
inline constexpr T PI (3.141592653589793238462643383279502884197169399375105820974944592307816406286L);

namespace knl {

  enum class coordsys_t : unsigned char
    {
     Cartesian = 0,
     Cylindrical,
     Spherical,
     LogSpherical,
     LogSphericalEV,
    };

  template < coordsys_t Coord = coordsys_t::Cartesian >
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

    template < class Vec_x, class Vec_v, typename T >
    static inline auto geodesic_move( Vec_x& x, const Vec_v& v, T dt ) {
      auto dx = v * dt;
      x += dx;
      return dx;
    }
  };

  template <>
  struct coord < coordsys_t::LogSpherical > {
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

    // template < typename T >
    // static inline Vec3<T> gdm_old( Vec3<T&>& x, Vec3<T&>& p, T dt ) {
    //   auto& logr  = std::get<0>(x);
    //   auto& theta = std::get<1>(x);
    //   auto& phi = std::get<2>(x);

    //   auto logr_o = logr;
    //   auto theta_o = theta;
    //   auto phi_o = phi;

    //   auto& p_r = std::get<0>(p);
    //   auto& p_t = std::get<1>(p);
    //   auto& p_ph = std::get<2>(p);

    //   auto r = std::exp(logr);
    //   auto sin_t = std::sin(theta);
    //   auto cos_t = std::cos(theta);
    //   auto sin_ph = std::sin(phi);
    //   auto cos_ph = std::cos(phi);

    //   // transform momentum to cartesian
    //   Vec3<T> p_cart ( (p_r * sin_t + p_t * cos_t) * cos_ph - p_ph * sin_ph,
    //                       (p_r * sin_t + p_t * cos_t) * sin_ph + p_ph * cos_ph,
    //                       p_r * cos_t - p_t * sin_t );
    //   // transform position to cartesian
    //   Vec3<T> x_cart ( r * sin_t * cos_ph,
    //                       r * sin_t * sin_ph,
    //                       r * cos_t );

    //   // update x_cart
    //   coord<coordsys_t::Cartesian>::geodesic_move( x_cart, p_cart, dt );

    //   // transform the new position back from cartesian
    //   logr = vecop::abs(x_cart);
    //   theta = std::acos(std::get<2>(x_cart) / logr);
    //   // TODO: Check correctness of this statement!
    //   phi = std::atan( std::get<1>(x_cart)/ std::get<0>(x_cart) ) + ( std::get<0>(x_cart) < 0 ) * PI<T>;
    //   logr = std::log(logr);

    //   // rebase momentum to the new position in the original coordinate system
    //   // These become vr, vtheta, vphi
    //   r = std::exp(logr);
    //   sin_t = std::sin(theta);
    //   cos_t = std::cos(theta);
    //   sin_ph = std::sin(phi);
    //   cos_ph = std::cos(phi);

    //   p_r = std::get<0>(p_cart) * sin_t * cos_ph + std::get<1>(p_cart) * sin_t * sin_ph + std::get<2>(p_cart) * cos_t;
    //   p_t = std::get<0>(p_cart) * cos_t * cos_ph + std::get<1>(p_cart) * cos_t * sin_ph - std::get<2>(p_cart) * sin_t;
    //   p_ph = -std::get<0>(p_cart) * sin_ph +  std::get<1>(p_cart) * cos_ph;

    //   return Vec3<T>( logr - logr_o, theta - theta_o, phi - phi_o);
    // }


    template < class Vec_x, class Vec_v, typename T >
    static inline auto geodesic_move( Vec_x& x, Vec_v& v, T dt ) {
      // TODOL: in implementing this, we assumed 1) no crossing through center and 2) no crossing through symmetry axes

      // dx is the return value and meanwhile it serves as a temporary
      auto dx = v * dt;
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

        auto[v_r, v_theta, v_phi] = v;

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
        phi -= 2 * PI<T> * std::floor( phi / 2 * PI<T> );
      }

      return dx;
    }
  };
}




#endif
