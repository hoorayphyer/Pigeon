#include "particle/forces.hpp"
#include "apt/numeric.hpp"

namespace particle {
  template < typename T, template < typename > class S >
  map<Force<T,S>> force_map;

  // TODO turn these into unique_ptr and use deepcopy
  template < typename T, template < typename > class S >
  void Force<T,S>::add ( force_t<T,S> force, T param ) {
    forces.push_back(force);
    params.push_back(param);
  }

  template < typename T, template < typename > class S >
  void Force<T,S>::Register( species sp ) const {
    force_map<T,S>[sp] = *this;
  }

  template < typename T, template < typename > class S >
  void Force<T,S>::Unregister( species sp ) {
    force_map<T,S>.erase(sp);
  }

  template < typename T, template < typename > class S >
  Force<T,S> Force<T,S>::Get( species sp ) noexcept {
    // if sp doesn't exist, newly created item in the maps will automatically have zero elements
    if ( force_map<T,S>.has(sp) ) return force_map<T,S>.at(sp);
    else return {};
  }
}

// TODO optimize use of intermediate variables
namespace particle::force {
  template < typename T, template < typename > class S, template < typename, template < typename > class > class Ptc_t >
  void lorentz( Ptc_t<T,S>& ptc, T dt, const apt::Vec<T,S<T>::Dim>& E, const apt::Vec<T,S<T>::Dim>& B, T q_over_m  ) {
    using Vec = apt::Vec<T,S<T>::Dim>;
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
