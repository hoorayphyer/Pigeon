#ifndef _PARTICLE_SCATTERING_HPP_
#define _PARTICLE_SCATTERING_HPP_

#include "particle/array.hpp"
#include "particle/properties.hpp"
#include "apt/vec.hpp"
#include "apt/numeric.hpp"
#include <optional>
#include <memory>
#include "random/rng.hpp"

namespace particle {
  template < typename T, template < typename > class S >
  using Ptc_t = typename array<T,S>::particle_type;

  namespace scat {
    template < typename T, template < typename > class S >
    using Eligible_t = bool (*) ( const Ptc_t<T,S>& );

    template < typename T, template < typename > class S >
    using Channel_t = std::optional<T> (*)( const Ptc_t<T,S>& ptc, const Properties& props, const apt::Vec<T,S<T>::Dim>& dp, T dt, const apt::Vec<T,S<T>::Dim>& B, util::Rng<T>& rng );
  }

  template < typename T, template < typename > class S >
  struct Scat {
    std::vector<scat::Eligible_t<T,S>> eligs;
    std::vector<scat::Channel_t<T,S>> channels;
    void (*impl) ( std::back_insert_iterator<array<T,S>> itr, Ptc_t<T,S>& ptc, T ) = nullptr;

    void Register( species sp ) const;
    static void Unregister( species sp );
    static Scat<T,S> Get ( species sp );
  };
}

// Channels
namespace particle::scat {
  template < typename T, template < typename > class S >
  struct CurvatureRadiation {
    using Vec = apt::Vec<T,S<T>::Dim>;

    static constexpr T calc_Rc ( const Ptc_t<T,S>& ptc, const Properties& props, const Vec& dp, T dt, const Vec& B ) noexcept {
      // TODOL uniform Rc is used here
      return 1.0;

      // qB / (\gamma m c) * dt < 2 \pi / 10
      bool is_gyration_resolved =
        std::sqrt( apt::sqabs(B) / ( (props.mass_x != 0) + apt::sqabs(ptc.p()) ) ) * std::abs(props.charge_x) * dt / props.mass_x < (2 * std::acos(-1) / 10.0);

      if ( is_gyration_resolved ) {
        // find momentum at half time step
        auto phalf = ptc.p() - dp * 0.5; // NOTE ptc.p() is already the updated p
        auto v = phalf / std::sqrt( (props.mass_x != 0) + apt::sqabs(phalf) );
        auto vv = apt::sqabs(v);
        Vec a = dp / dt; // a is for now force, will be converted to dv/dt
        // convert a to dv/dt
        a = ( a - v * apt::dot(v,a) ) * std::sqrt( 1.0 - vv );
        auto va = apt::dot(v,a); // get the real v dot a
        return vv / std::max( std::sqrt( apt::sqabs(a) - va * va / vv ), 1e-6 ); // in case denominator becomes zero
      } else {
        // Dipolar radius of curvature in LogSpherical
        const auto& theta = ptc.q()[1];
        auto tmp = 2.5 + 1.5 * std::cos( 2 * theta );
        return std::exp(ptc.q()[0]) * tmp * std::sqrt(tmp) / ( ( tmp + 2.0 ) * std::sin(theta) );
      }
    }

    static T K_thr;
    static T gamma_off;
    static T emission_rate;
    // TODO move calc_Rc out too
    // static T (*calc_Rc) ( const Ptc_t<T,S>& ptc, const Properties& props, const Vec& dp, T dt, const Vec& B );
    static T (*sample_E_ph)();

    // return sampled energy if any
    static constexpr std::optional<T> test ( const Ptc_t<T,S>& ptc, const Properties& props, const Vec& dp, T dt,
                                             const Vec& B, util::Rng<T>& rng ) noexcept {
      T Rc = calc_Rc( ptc, props, dp, dt, B );
      T gamma = std::sqrt( (props.mass_x != 0) + apt::sqabs(ptc.p()) );

      if (gamma > gamma_off && gamma > K_thr *  std::cbrt(Rc) && rng.uniform() < emission_rate * dt * gamma  / Rc ) {
        return { std::min<T>( sample_E_ph(), gamma - 1.0 ) };
      } else return {};
    }
  };

  template < typename T, template < typename > class S >
  struct MagneticConvert {
    static T B_thr;
    static T mfp;

    static constexpr std::optional<T> test ( const Ptc_t<T,S>& photon, const Properties&, const apt::Vec<T,S<T>::Dim>& , T dt,
                                         const apt::Vec<T,S<T>::Dim>& B, util::Rng<T>& rng ) noexcept {
      return ( apt::sqabs(B) > B_thr * B_thr ) && ( rng.uniform() < dt / mfp )  ? std::optional<T>(1.0) : std::nullopt; // NOTE return std::optional(1.0) so as to be treated true in boolean conversion
    }
  };

  template < typename T, template < typename > class S >
  struct TwoPhotonCollide {
    // TODO double check this implementation, it is not equivalent because the original one has some sort of gaussian in it. Use Monte Carlo
    // inline T f_x ( T x ) {
    //   // distribution of x*exp(-x^2/2), which peaks at x = 1.
    //   return std::sqrt( -2.0 * std::log(x) );
    // }
    static T mfp;

    static constexpr std::optional<T> test ( const Ptc_t<T,S>& , const Properties&, const apt::Vec<T,S<T>::Dim>& , T dt,
                                                   const apt::Vec<T,S<T>::Dim>& , util::Rng<T>& rng ) noexcept {
      return ( rng.uniform() < dt / mfp )  ? std::optional<T>(1.0) : std::nullopt;
    }
  };
}

// Impls
namespace particle::scat {
  template < bool Instant, typename T, template < typename > class S >
  void RadiationFromCharges ( std::back_insert_iterator<array<T,S>> itr, Ptc_t<T,S>& ptc, T E_ph );

  template < typename T, template < typename > class S >
  void PhotonPairProduction ( std::back_insert_iterator<array<T,S>> itr, Ptc_t<T,S>& photon, T );
}

#endif
