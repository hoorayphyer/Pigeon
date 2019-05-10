#include "particle/scattering.hpp"
#include "particle/particle.hpp"
#include "apt/numeric.hpp"

namespace particle::scat {

  template < typename T, template < typename > class PtcSpecs >
  using real_Ptc = Particle<T,PtcSpecs>;

  template < bool Instant, typename T, template < typename > class PtcSpecs >
  void RadiationFromCharges<Instant,T,PtcSpecs>::impl ( std::back_insert_iterator<array<T,PtcSpecs>> itr, Ptc_t<T,PtcSpecs>& ptc, T E_ph ) {
    {
      auto gamma_ptc = std::sqrt( 1.0 + apt::sqabs(ptc.p()) );
      // primary particle loses energy to gamma rays. ptc.p *= |pf| / |pi|
      ptc.p() *= std::sqrt( ( 1.0 - E_ph / ( gamma_ptc - 1.0 ) ) * ( 1.0 - E_ph / ( gamma_ptc + 1.0 ) ) );
    }

    if constexpr ( Instant ) {
        // recycle Rc for gamma_sec
        E_ph /= 2.0;
        // append electron and positron
        auto ptc_sec = real_Ptc<T,PtcSpecs> ( ptc.q(),
                                              ptc.p() * ( std::sqrt( E_ph * E_ph - 1.0 ) / apt::sqabs(ptc.p()) ),
                                              flag::secondary );
        ptc_sec.set(species::electron);
        *(itr++) = ptc_sec;
        ptc_sec.set(species::positron);
        *(itr++) = std::move(ptc_sec);

      } else {
      *(itr++) = real_Ptc<T,PtcSpecs> ( ptc.q(),
                                        ptc.p() * ( E_ph / apt::abs(ptc.p()) ),
                                        species::photon );
    }
  }

  template < typename T, template < typename > class PtcSpecs >
  void PhotonPairProduction<T, PtcSpecs>::impl( std::back_insert_iterator<array<T,PtcSpecs>> itr, Ptc_t<T,PtcSpecs>& photon, T ) {
    auto ptc_sec = real_Ptc<T,PtcSpecs> ( photon.q(),
                                          photon.p() * std::sqrt( 0.25 - 1.0 / apt::sqabs(photon.p()) ),
                                          flag::secondary );
    ptc_sec.set(species::electron);
    *(itr++) = ptc_sec;
    ptc_sec.set(species::positron);
    *(itr++) = std::move(ptc_sec);

    // void this photon
    photon.set(flag::empty);

  }
}


#include "particle/map.hpp"

namespace particle {
  template < typename T, template < typename > class PtcSpecs >
  map<std::unique_ptr<scat::Scat<T,PtcSpecs>>> scat_map;

  template < typename T, template < typename > class PtcSpecs >
  void ScatGen<T,PtcSpecs>::Register( species sp, const scat::Scat<T,PtcSpecs>& scat ) {
    scat_map<T,PtcSpecs>[sp].reset( new scat::Scat<T,PtcSpecs>(scat) ); // NOTE use copy constructor
  }

  template < typename T, template < typename > class PtcSpecs >
  void ScatGen<T,PtcSpecs>::Unregister( species sp ) {
    scat_map<T,PtcSpecs>.erase(sp);
  }

  template < typename T, template < typename > class PtcSpecs >
  scat::Scat<T,PtcSpecs>* ScatGen<T,PtcSpecs>::operator() ( species sp ) {
    return scat_map<T,PtcSpecs>.has(sp) ? scat_map<T,PtcSpecs>.at(sp).get() : nullptr;
  }
}
