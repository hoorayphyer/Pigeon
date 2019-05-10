#include "particle/forces.hpp"
#include "apt/numeric.hpp"

// TODO optimize use of intermediate variables
namespace particle::force {
  template < typename T, template < typename > class PtcSpecs, template < typename, template < typename > class > class Ptc_t >
  void lorentz( Ptc_t<T,PtcSpecs>& ptc, T dt, const Vec<T,PtcSpecs>& E, const Vec<T,PtcSpecs>& B, T q_over_m  ) {
    using Vec = Vec<T,PtcSpecs>;
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
}
