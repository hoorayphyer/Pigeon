#include "particle/scattering_impl.hpp"
#include "pic.hpp"

namespace particle {
  using namespace pic;

  template <>
  real_t scat::CurvatureRadiation<real_t,Specs>::K_thr{};
  template <>
  real_t scat::CurvatureRadiation<real_t,Specs>::gamma_off{};
  template <>
  real_t scat::CurvatureRadiation<real_t,Specs>::emission_rate{};
  template <>
  real_t (*scat::CurvatureRadiation<real_t,Specs>::sample_E_ph) () = nullptr;


  template <>
  real_t scat::MagneticConvert<real_t,Specs>::B_thr {};
  template <>
  real_t scat::MagneticConvert<real_t,Specs>::mfp {};

  template <>
  real_t scat::TwoPhotonCollide<real_t,Specs>::mfp {};


  namespace scat {
    template void RadiationFromCharges<true>(std::back_insert_iterator<array<real_t,Specs>>, Ptc_t<real_t,Specs>&, real_t);
    template void RadiationFromCharges<false>(std::back_insert_iterator<array<real_t,Specs>>, Ptc_t<real_t,Specs>&, real_t);
    template void PhotonPairProduction(std::back_insert_iterator<array<real_t,Specs>>, Ptc_t<real_t,Specs>&, real_t);
  }
  template class ScatGen<real_t, Specs>;
}
