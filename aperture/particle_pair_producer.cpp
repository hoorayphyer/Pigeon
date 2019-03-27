#include "particle/pair_producer.cpp"
#include "traits.hpp"
namespace particle {
  using namespace traits;
  template struct LeptonProduces<true, real_t, DPtc, ptc_state_t>;
  template struct LeptonProduces<false, real_t, DPtc, ptc_state_t>;
  template struct PhotonProduces<real_t, DPtc, ptc_state_t>;
}
