#include "particle/scattering.hpp"
#include "particle/particle.hpp"
#include "apt/numeric.hpp"

namespace particle::scat {

  template < class PtcArr >
  using real_Ptc = Particle<typename PtcArr::value_type, PtcArr::DPtc, typename PtcArr::state_type>;

  template < bool Instant, class PtcArr >
  void RadiationFromCharges<Instant,PtcArr>::impl ( std::back_insert_iterator<PtcArr> itr, ts::Ptc<PtcArr>& ptc, ts::T<PtcArr> E_ph ) {
    {
      auto gamma_ptc = std::sqrt( 1.0 + apt::sqabs(ptc.p()) );
      // primary particle loses energy to gamma rays. ptc.p *= |pf| / |pi|
      ptc.p() *= std::sqrt( ( 1.0 - E_ph / ( gamma_ptc - 1.0 ) ) * ( 1.0 - E_ph / ( gamma_ptc + 1.0 ) ) );
    }

    if constexpr ( Instant ) {
        // recycle Rc for gamma_sec
        E_ph /= 2.0;
        // append electron and positron
        auto ptc_sec = real_Ptc<PtcArr> ( ptc.q(),
                                          ptc.p() * ( std::sqrt( E_ph * E_ph - 1.0 ) / apt::sqabs(ptc.p()) ),
                                          flag::secondary );
        ptc_sec.set(species::electron);
        *(itr++) = ptc_sec;
        ptc_sec.set(species::positron);
        *(itr++) = std::move(ptc_sec);

      } else {
      *(itr++) = real_Ptc<PtcArr> ( ptc.q(),
                                    ptc.p() * ( E_ph / apt::abs(ptc.p()) ),
                                    species::photon );
    }
  }

  template < class PtcArr >
  void PhotonPairProduction<PtcArr>::impl( std::back_insert_iterator<PtcArr> itr, ts::Ptc<PtcArr>& photon, ts::T<PtcArr> ) {
    auto ptc_sec = real_Ptc<PtcArr> ( photon.q(),
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
  template < class PtcArr >
  map<std::unique_ptr<scat::Scat<PtcArr>>> scat_map;

  template < class PtcArr >
  void ScatGen<PtcArr>::Register( species sp, const scat::Scat<PtcArr>& scat ) {
    scat_map<PtcArr>[sp].reset( new scat::Scat<PtcArr>(scat) ); // NOTE use copy constructor
  }

  template < class PtcArr >
  void ScatGen<PtcArr>::Unregister( species sp ) {
    scat_map<PtcArr>.erase(sp);
  }

  template < class PtcArr >
  scat::Scat<PtcArr>* ScatGen<PtcArr>::operator() ( species sp ) {
    return scat_map<PtcArr>.has(sp) ? scat_map<PtcArr>.at(sp).get() : nullptr;
  }
}
