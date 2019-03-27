#ifndef _PAIR_PRODUCER_HPP_
#define _PAIR_PRODUCER_HPP_

#include "particle/array.hpp"
#include "apt/vec.hpp"
#include <optional>

namespace particle {
  template < typename T, int DPtc, typename state_t >
  struct EmittedEnergy {
    std::optional<T> operator() ( const vParticle<T, DPtc, state_t>& ptc, const apt::Vec<T,DPtc>& dp ) = 0;
  };

}

namespace particle {
  template < bool Instant, typename T, int DPtc, typename state_t >
  struct LeptonProduces {
    void operator() ( std::back_insert_iterator<array<T,DPtc,state_t>> itr, vParticle<T,DPtc,state_t>& ptc, T emitted_energy );
  };

  template < typename T, int DPtc, typename state_t >
  struct PhotonProduces {
    void operator() ( std::back_insert_iterator<array<T,DPtc,state_t>> itr, vParticle<T,DPtc,state_t>& photon );
  };
}

#endif
