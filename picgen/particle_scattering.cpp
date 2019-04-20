#include "particle/scattering.cpp"
#include "pic.hpp"

namespace particle {
  using namespace pic;

  namespace scat {
    template class RadiationFromCharges<true, real_t, Specs>;
    template class RadiationFromCharges<false, real_t, Specs>;
    template class PhotonPairProduction<real_t, Specs>;
  }
  template class ScatGen<real_t, Specs>;
}
