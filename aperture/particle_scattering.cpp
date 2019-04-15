#include "particle/scattering.cpp"
#include "traits.hpp"
namespace particle {
  using namespace traits;
  using PtcArr = array<real_t,  DPtc, ptc_state_t>;

  namespace scat {
    template class RadiationFromCharges<true, PtcArr>;
    template class RadiationFromCharges<false, PtcArr>;
    template class PhotonPairProduction<PtcArr>;
  }
  template class ScatGen<PtcArr>;
}
