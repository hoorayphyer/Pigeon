#include "particle/scattering.hpp"
#include "particle/particle.hpp"
#include "apt/numeric.hpp"

namespace particle::scat {
  template < bool Instant, typename T, template < typename > class S >
  void RadiationFromCharges ( std::back_insert_iterator<array<T,S>> itr, Ptc_t<T,S>& ptc, T E_ph ) {
    {
      auto gamma_ptc = std::sqrt( 1.0 + apt::sqabs(ptc.p()) );
      // primary particle loses energy to gamma rays. ptc.p *= |pf| / |pi|
      ptc.p() *= std::sqrt( ( 1.0 - E_ph / ( gamma_ptc - 1.0 ) ) * ( 1.0 - E_ph / ( gamma_ptc + 1.0 ) ) );
    }

    if constexpr ( Instant ) {
        // recycle Rc for gamma_sec
        E_ph /= 2.0;
        // append electron and positron
        Particle<T,S> ptc_sec ( ptc.q(),
                                ptc.p() * ( std::sqrt( E_ph * E_ph - 1.0 ) / apt::sqabs(ptc.p()) ),
                                flag::secondary,
                                ptc.template get<birthplace>() );
        ptc_sec.set(species::electron);
        *(itr++) = ptc_sec;
        ptc_sec.set(species::positron);
        *(itr++) = std::move(ptc_sec);

      } else {
      *(itr++) = Particle<T,S> ( ptc.q(),
                                 ptc.p() * ( E_ph / apt::abs(ptc.p()) ),
                                 species::photon,
                                 ptc.template get<birthplace>()
                                 );
    }
  }

  template < typename T, template < typename > class S >
  void PhotonPairProduction ( std::back_insert_iterator<array<T,S>> itr, Ptc_t<T,S>& photon, T ) {
    Particle<T,S> ptc_sec ( photon.q(),
                            photon.p() * std::sqrt( 0.25 - 1.0 / apt::sqabs(photon.p()) ),
                            flag::secondary,
                            photon.template get<birthplace>()
                            );
    ptc_sec.set(species::electron);
    *(itr++) = ptc_sec;
    ptc_sec.set(species::positron);
    *(itr++) = std::move(ptc_sec);

    // void this photon
    photon.reset(flag::exist);
  }
}

#include "particle/map.hpp"

namespace particle {
  template < typename T, template < typename > class S >
  map<scat::Scat<T,S>> scat_map;

  template < typename T, template < typename > class S >
  void ScatGen<T,S>::Register( species sp, const scat::Scat<T,S>& scat ) {
    scat_map<T,S>[sp] = scat;
  }

  template < typename T, template < typename > class S >
  void ScatGen<T,S>::Unregister( species sp ) {
    scat_map<T,S>.erase(sp);
  }

  template < typename T, template < typename > class S >
  scat::Scat<T,S>* ScatGen<T,S>::operator() ( species sp ) {
    return scat_map<T,S>.has(sp) ? &(scat_map<T,S>.at(sp)) : nullptr;
  }
}
