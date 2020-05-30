#include "particle/scattering_impl.hpp"
#include "pic.hpp"

namespace particle {
  using namespace pic;

  template <>
  real_t scat::CurvatureRadiation<real_t,Specs>::gamma_fd{};
  template <>
  real_t scat::CurvatureRadiation<real_t,Specs>::gamma_off{};
  template <>
  real_t scat::CurvatureRadiation<real_t,Specs>::Ndot_fd{};
  template <>
  real_t scat::CurvatureRadiation<real_t,Specs>::E_ph{};
  // template <>
  // real_t (*scat::CurvatureRadiation<real_t,Specs>::sample_E_ph) () = nullptr;


  template <>
  real_t scat::MagneticConvert<real_t,Specs>::B_thr {};
  template <>
  real_t scat::MagneticConvert<real_t,Specs>::mfp {};

  template <>
  real_t scat::TwoPhotonCollide<real_t,Specs>::mfp {};


  namespace scat {
    template void RadiationFromCharges<true>(back_insert_iterator_t<real_t,Specs>, Ptc_t<real_t,Specs>&, real_t, flagger_t);
    template void RadiationFromCharges<false>(back_insert_iterator_t<real_t,Specs>, Ptc_t<real_t,Specs>&, real_t, flagger_t);
    template void PhotonPairProduction(back_insert_iterator_t<real_t,Specs>, Ptc_t<real_t,Specs>&, real_t, flagger_t);
  }
  template class Scat<real_t, Specs>;
}
