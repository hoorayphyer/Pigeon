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
  void lorentz( Ptc_t<T,S>& ptc, T dt, const apt::Vec<T,S<T>::Dim>& E, const apt::Vec<T,S<T>::Dim>& B, T q_times_w_gyro_unitB_over_m  ) {
    using Vec = apt::Vec<T,S<T>::Dim>;
    // lambda = 0.5 * w_gyro_unitB * dt * (charge_x unit_q) / (mass_x unit_m) NOTE this is actually rescaling Lorentz force
    dt *= 0.5 * q_times_w_gyro_unitB_over_m; // repurpose dt for lambda

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

  // NOTE this is the Taylor expansion of (1/x) * ln( ( x+b+sqrt(1+2bx+x^2) ) / ( 1+b ) )
  // NOTE this requires e E dt /m c << 1
  template < typename T >
  constexpr T F( T x, T b ) noexcept {
    return 1 - 0.5 * b * x + (1.0 / 6.0) * ( 3 * b * b - 1 ) * x * x;
  }

  template < typename T, template < typename > class S, template < typename, template < typename > class > class Ptc_t >
  void lorentz_exact( Ptc_t<T,S>& ptc, T dt, const apt::Vec<T, S<T>::Dim>& E, const apt::Vec<T, S<T>::Dim>& B, T q_times_w_gyro_unitB_over_m  ) {
    using Vec = apt::Vec<T,S<T>::Dim>;
    using vVec = apt::vVec<T,S<T>::Dim>;

    // get drift velocity
    T xp = apt::dot(E,B);
    T xE = apt::sqabs(B) - apt::sqabs(E); // calculate B^2 - E^2
    T xB = ( std::sqrt(xE * xE + 4 * xp * xp) + xE ) / 2.0; // calculate B'^2

    Vec beta_d = apt::cross(E,B) / ( apt::sqabs(E) + xB );
    T gamma_d = 1.0 / std::sqrt( 1 - apt::sqabs(beta_d) );

    T gamma {}, gamma_co{};
    { /// ------- boost particle momentum to comoving frame
      gamma = sqrt( 1 + apt::sqabs(ptc.p()) );
      gamma_co = gamma_d * ( gamma - apt::dot(beta_d, ptc.p()) );
      ptc.p() -= beta_d * ( gamma_d * ( gamma + gamma_co ) / ( 1.0 + gamma_d ) );
    }

    { /// ------- update in comoving frame
      Vec z;
      if ( xE > 0 ) { // test if B' > E'
        z = B - apt::cross(beta_d,E);
        z /= apt::abs(z);
        xE = ( ( xp > 0 ) - ( xp < 0 ) ) * sqrt( xB - xE ) * q_times_w_gyro_unitB_over_m * dt;
        xB = std::sqrt(xB) * q_times_w_gyro_unitB_over_m * dt;
      } else {
        z = E + apt::cross(beta_d,B);
        z /= apt::abs(z);
        xE = sqrt( xB - xE ) * q_times_w_gyro_unitB_over_m * dt;
        xB = ( ( xp > 0 ) - ( xp < 0 ) ) * std::sqrt(xB) * q_times_w_gyro_unitB_over_m * dt;
      }

      xp = apt::dot(z,ptc.p());

      gamma_co = std::sqrt( gamma_co * gamma_co + 2 * xp * 0.5 * xE + 0.25 * xE * xE );
      xp += 0.5 * xE; // 1st half update
      // update p_z
      ptc.p() += xE * z;

      xE *= 0.5;
      xB *= F(xE/gamma_co,xp/gamma_co) / gamma_co;

      xp += 0.5 * 2.0 * xE; // 2nd half update

      { // rotate p around B in the comoving frame
        ptc.p() = xp * z + std::cos(xB) * ( ptc.p() - xp * z);
        // by now, xp and xE are free
        xp = z[1] * ptc.p()[2] - z[2] * ptc.p()[1]; // ( z x p )_0
        xE = z[2] * ptc.p()[0] - z[0] * ptc.p()[2]; // ( z x p )_1
        z[2] = z[0] * ptc.p()[1] - z[1] * ptc.p()[0]; // ( z x p )_2
        z[0] = std::tan(xB);
        ptc.p()[0] += z[0] * xp;
        ptc.p()[1] += z[0] * xE;
        ptc.p()[2] += z[0] * z[2];
      }

      // TODO udate delta q
    }

    { /// ------- boost particle momentum back to lab
      gamma_co = sqrt( 1 + apt::sqabs(ptc.p()) );
      gamma = gamma_d * ( gamma_co + apt::dot(beta_d, ptc.p()) );
      ptc.p() += beta_d * ( gamma_d * ( gamma + gamma_co ) / ( 1.0 + gamma_d ) );
    }

  }
}
