#include "particle/pusher.hpp"
#include "particle/properties.hpp"
#include "apt/numeric.hpp"
#include "kernel/coordinate.hpp"

namespace particle :: force {
  auto landau0 =
    [](const auto& p, const auto& E, const auto& B) {
      auto EB2 = apt::dot(E,B);
      EB2 = EB2 * EB2;
      auto B2_E2 = apt::sqabs(B) - apt::sqabs(E);
      // calculate E'^2
      auto Ep2 = 2 * EB2 / ( std::sqrt(B2_E2 * B2_E2 + 4 * EB2) + B2_E2 );
      auto beta_ExB = apt::cross(E,B) / ( apt::sqabs(B) + Ep2);
      // find B' modulo gamma_ExB
      auto Bp = B - apt::cross( beta_ExB, E);
      // obtain the momentum with perpendicular components damped
      auto p_new = Bp * ( apt::dot( p, Bp ) / apt::sqabs(Bp) );
      p_new += beta_ExB * std::sqrt( ( 1.0 + apt::sqabs(p_new) ) / ( 1.0 - apt::sqabs(beta_ExB) ) );
      return p_new - p;
    };

  // lambda = dt / mass_x * e/m NOTE this is actually rescaling Lorentz force
  template <typename E_p, typename E_E, typename E_B, typename T = typename E_p::value_type, int N = E_p::size>
  apt::Vec<T,N> lorentz ( T lambda, const apt::VecExpression<E_p>& p, const apt::VecExpression<E_E>& E, const apt::VecExpression<E_B>& B ) noexcept {
    using Vec = apt::Vec<T,N>;
    // TODO optimize use of intermediate variables
    lambda /= 2.0;

    Vec u_halfstep = p + E * lambda + apt::cross(p, B) * ( lambda / std::sqrt( 1.0 + apt::sqabs(p) ) );
    Vec upr = u_halfstep + E * lambda;
    Vec tau = B * lambda;
    // store some repeatedly used intermediate results
    auto tt = apt::sqabs(tau);
    auto ut = apt::dot(upr, tau);

    auto sigma = 1.0 + apt::sqabs(upr) - tt;
    auto inv_gamma2 =  2.0 / ( sigma + std::sqrt( sigma * sigma + 4.0 * ( tt + ut * ut ) ) ); // inv_gamma2 means ( 1 / gamma^(i+1) ) ^2
    auto s = 1.0 / ( 1.0 + inv_gamma2 * tt );
    Vec p_vay = ( upr + tau * ( ut * inv_gamma2 ) + apt::cross(upr, tau) * std::sqrt(inv_gamma2) ) * s;
    return Vec(p_vay - p);
  };

}

// TODO move check of forces on_off to somewhere else
namespace particle {
  template < species sp, typename dp_t, typename Ptc, typename Field, typename T >
  dp_t update_p( PtcExpression<Ptc>& ptc, const T& dt,
                 const apt::VecExpression<Field>& E, const apt::VecExpression<Field>& B ) {
    dp_t dp;

    dp += force::lorentz( dt / mass_x<sp>, ptc.p(), E, B );
    // TODO pane
    // Apply Lorentz force
    // if ( _pane.lorentz_On  ) {
    //   dp += force::lorentz( dt / mass_x<sp>, ptc.p, E, B );
    // }

    // if ( _pane.gravity_On )
    //   dp += _pane.gravity( ptc.q ) * dt; // TODO gravity interface

    // // FIXME add control in dashboard for rad_cooling
    // // if ( is_radiative<sp> ) {
    // //   rad_cooling(p, BVector, EVector, 2e-9, dt);
    // // }

    // // when B is strong enough, damp the perpendicular component of momentum
    // // FIXME ions should also be affected right?
    // // TODO should this go before p += dp?
    // if ( _pane.landau0_On && is_radiative<sp> ) {
    //   if ( apt::sqabs(B) > _pane.B_landau0 * _pane.B_landau0 ) {
    //     dp += force::landau0( ptc.p, E, B );
    //   }
    // }

    ptc.p() += dp;

    return dp;
  }


  template < species sp, knl::coordsys CS, typename dq_t, typename Ptc, typename T >
  dq_t update_q( PtcExpression<Ptc>& ptc, const T& dt ) {
    auto gamma = std::sqrt( is_massive<sp> + apt::sqabs(ptc.p()) );

    if constexpr ( CS == knl::coordsys::Cartesian ) {
      // a small optimization for Cartesian
      return knl::coord<CS>::geodesic_move( ptc.q(), ptc.p(), dt / gamma );
    } else {
      dq_t dq = knl::coord<CS>::geodesic_move( ptc.q(), (ptc.p() /= gamma), dt );
      ptc.p() *= gamma;
      return dq;
    }

  }


}
