#include "particle_module/pusher.hpp"
#include "core/particle.hpp"
#include "core/coordinate.hpp"

namespace vn = apt::numerical;

namespace particle :: force {
  auto landau0 =
    [](const auto& p, const auto& E, const auto& B) {
      Real EB2 = vn::E.dot(B);
      EB2 = EB2 * EB2;
      Real B2_E2 = vn::abs_sq(B) - vn::abs_sq(E);
      // calculate E'^2
      Real Ep2 = 2 * EB2 / ( std::sqrt(B2_E2 * B2_E2 + 4 * EB2) + B2_E2 );
      auto beta_ExB = vn::cross(E,B) / ( vn::abs_sq(B) + Ep2);
      // find B' modulo gamma_ExB
      auto Bp = B - vn::cross( beta_ExB, E);
      // obtain the momentum with perpendicular components damped
      auto p_new = Bp * ( vn::dot( p, Bp ) / vn::abs_sq(Bp) );
      p_new += beta_ExB * std::sqrt( ( 1.0 + vn::abs_sq(p_new) ) / ( 1.0 - vn::abs_sq(beta_ExB) ) );
      return p_new - p;
    }

  auto lorentz = // lambda = dt / mass * e/m NOTE this is actually rescaling Lorentz force
    []( Real lambda, const auto& p, const auto& E, const auto& B ) noexcept {
      // TODO optimize use of intermediate variables
      lambda /= 2.0;

      auto u_halfstep = p + E * lambda + vn::cross(p, B) * ( lambda / std::sqrt( 1.0 + vn::abs_sq(p) ) );
      auto upr = u_halfstep + E * lambda;
      auto tau = B * lambda;
      // store some repeatedly used intermediate results
      Real tt = vn::abs_sq(tau);
      Real ut = vn::dot(upr, tau);

      Real sigma = 1.0 + vn::abs_sq(upr) - tt;
      Real inv_gamma2 =  2.0 / ( sigma + std::sqrt( sigma * sigma + 4.0 * ( tt + ut * ut ) ) ); // inv_gamma2 means ( 1 / gamma^(i+1) ) ^2
      Real s = 1.0 / ( 1.0 + inv_gamma2 * tt );
      auto p_vay = ( upr + tau * ( ut * inv_gamma2 ) + vn::cross(upr, tau) * std::sqrt(inv_gamma2) ) * s;
      return p_vay - p;
    };

}

// TODO move check of forces on_off to somewhere else
namespace particle {
  template < typename Tvt, std::size_t DPtc, std::size_t DField,
             typename Trl = apt::remove_cvref_t<Tvt> >
  Vec<Trl,DPtc> update_p( Particle<Tvt,DPtc>& ptc, const Species& sp, Trl dt,
                          const Vec<Trl, DField>& E, const Vec<Trl, DField>& B ) {
    Vec<Trl, DPtc> dp;

    // Apply Lorentz force
    if ( _pane.lorentz_On  ) {
      dp += force::lorentz( dt / sp.mass, ptc.p, E, B );
    }

    if ( _pane.gravity_On )
      dp += _pane.gravity( ptc.q ) * dt; // TODO gravity interface

    // TODO put Rc in a separate container
    // // FIXME don't use a uniform number for Rc
    // // ptc.Rc = CalculateRc( dt, p, dp );
    // ptc.Rc = 1.0;

    // FIXME add control in dashboard for rad_cooling
    // if ( sp.is_radiative ) {
    //   rad_cooling(p, BVector, EVector, 2e-9, dt);
    // }

    // when B is strong enough, damp the perpendicular component of momentum
    // FIXME ions should also be affected right?
    // TODO should this go before p += dp
    if ( _pane.landau0_On && sp.is_radiative ) {
      if ( vn::abs_sq(B) > _pane.B_landau0 * _pane.B_landau0 ) {
        dp += force::landau0( ptc.p, E, B );
      }
    }

    ptc.p += dp;

    return dp;
  }


  template < CoordSys CS, typename Tvt, std::size_t DPtc,
             typename Trl = apt::remove_cvref_t<Tvt> >
  Vec<Trl,DPtc> update_q( Particle<Tvt, DPtc>& ptc, const Species& sp, Trl dt ) {
    Real gamma = std::sqrt( (sp.mass > 0) + vn::abs_sq(ptc.p) );

    if constexpr ( CS == CoordSys::Cartesian ) {
        return coord<CS>::geodesic_move( ptc.q, ptc.p, dt / gamma );
      } else {
      auto dq = coord<CS>::geodesic_move( ptc.q, (ptc.p /= gamma), dt );
      ptc.p *= gamma;
      return dq;
    }

  }


}
