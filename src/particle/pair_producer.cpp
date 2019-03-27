#include "particle/pair_producer.hpp"
#include "particle/particle.hpp"
#include "apt/numeric.hpp"

namespace particle {
  template < bool Instant, typename T, int DPtc, typename state_t >
  void LeptonProduces<Instant,T,DPtc,state_t>::operator()( std::back_insert_iterator<array<T,DPtc,state_t>> itr, vParticle<T,DPtc,state_t>& ptc, T E_ph ) {
    {
      T gamma_ptc = std::sqrt( 1 + apt::sqabs(ptc.p()) );
      // primary particle loses energy to gamma rays. ptc.p *= |pf| / |pi|
      ptc.p() *= std::sqrt( ( 1.0 - E_ph / ( gamma_ptc - 1.0 ) ) * ( 1.0 - E_ph / ( gamma_ptc + 1.0 ) ) );
    }

    if constexpr ( Instant ) {
        // recycle Rc for gamma_sec
        E_ph /= 2.0;
        // append electron and positron
        auto ptc_sec = Particle<T, DPtc, state_t> ( ptc.q(),
                                                    ptc.p() * ( std::sqrt( E_ph * E_ph - 1.0 ) / apt::sqabs(ptc.p()) ),
                                                    flag::secondary );
        ptc_sec.set(species::electron);
        *(itr++) = ptc_sec;
        ptc_sec.set(species::positron);
        *(itr++) = std::move(ptc_sec);

      } else {
      *(itr++) = Particle<T, DPtc, state_t > ( ptc.q(),
                                               ptc.p() * ( E_ph / apt::abs(ptc.p()) ),
                                               species::photon );
    }
  }

  template < typename T, int DPtc, typename state_t >
  void PhotonProduces<T,DPtc,state_t>::operator()( std::back_insert_iterator<array<T,DPtc,state_t>> itr, vParticle<T,DPtc,state_t>& photon ) {
    auto ptc_sec = Particle<T, DPtc, state_t> ( photon.q(),
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
