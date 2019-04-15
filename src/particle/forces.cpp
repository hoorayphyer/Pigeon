#include "particle/forces.hpp"
#include "apt/numeric.hpp"

// TODO optimize use of intermediate variables
namespace particle::force {
  template < typename Ptc >
  void lorentz( Ptc& ptc, ts::Real<Ptc> dt, const ts::Vec<Ptc>& E, const ts::Vec<Ptc>& B, ts::Real<Ptc> q_over_m  ) {
    using Vec = ts::Vec<Ptc>;
    // lambda = 0.5 * dt * (charge_x unit_q) / (mass_x unit_m) NOTE this is actually rescaling Lorentz force
    dt *= 0.5 * q_over_m; // repurpose dt for lambda

    const auto& p = ptc.p();
    Vec u_halfstep = p + E * dt + apt::cross(p, B) * ( dt / std::sqrt( 1.0 + apt::sqabs(p) ) );
    Vec upr = u_halfstep + E * dt;
    Vec tau = B * dt;
    // store some repeatedly used intermediate results
    auto tt = apt::sqabs(tau);
    auto ut = apt::dot(upr, tau);

    auto sigma = 1.0 + apt::sqabs(upr) - tt;
    auto inv_gamma2 =  2.0 / ( sigma + std::sqrt( sigma * sigma + 4.0 * ( tt + ut * ut ) ) ); // inv_gamma2 means ( 1 / gamma^(i+1) ) ^2
    auto s = 1.0 / ( 1.0 + inv_gamma2 * tt );
    ptc.p() = ( upr + tau * ( ut * inv_gamma2 ) + apt::cross(upr, tau) * std::sqrt(inv_gamma2) ) * s;
  }

  // when B is strong enough, damp the perpendicular component of momentum
  template < typename Ptc >
  void landau0( Ptc& ptc, ts::Real<Ptc> dt, const ts::Vec<Ptc>& E, const ts::Vec<Ptc>& B, ts::Real<Ptc> B2_thr  ) {
    if ( apt::sqabs(B) < B2_thr ) return;

    auto EB2 = apt::dot(E,B);
    EB2 = EB2 * EB2;
    auto B2_E2 = apt::sqabs(B) - apt::sqabs(E);
    // calculate E'^2
    auto Ep2 = 2 * EB2 / ( std::sqrt(B2_E2 * B2_E2 + 4 * EB2) + B2_E2 );
    auto beta_ExB = apt::cross(E,B) / ( apt::sqabs(B) + Ep2);
    // find B' modulo gamma_ExB
    auto Bp = B - apt::cross( beta_ExB, E);
    // obtain the momentum with perpendicular components damped
    ptc.p() = Bp * ( apt::dot( ptc.p(), Bp ) / apt::sqabs(Bp) );
    ptc.p() += beta_ExB * std::sqrt( ( 1.0 + apt::sqabs(ptc.p()) ) / ( 1.0 - apt::sqabs(beta_ExB) ) );
  };
}
