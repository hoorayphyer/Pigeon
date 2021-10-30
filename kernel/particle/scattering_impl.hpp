#include "apt/numeric.hpp"
#include "particle/particle.hpp"
#include "particle/scattering.hpp"

namespace particle::scat {
template <bool Instant, typename T, template <typename> class S>
void RadiationFromCharges(back_insert_iterator_t<T, S> itr, Ptc_t<T, S>& ptc,
                          T E_ph, flagger_t f) {
  {
    auto gamma_ptc = std::sqrt(1.0 + apt::sqabs(ptc.p()));
    // primary particle loses energy to gamma rays. ptc.p *= |pf| / |pi|
    ptc.p() *= std::sqrt((1.0 - E_ph / (gamma_ptc - 1.0)) *
                         (1.0 - E_ph / (gamma_ptc + 1.0)));
  }

  if constexpr (Instant) {
    // recycle Rc for gamma_sec
    E_ph /= 2.0;
    // append electron and positron
    Particle<T, S> ptc_sec(
        ptc.q(), ptc.p() * (std::sqrt(E_ph * E_ph - 1.0) / apt::sqabs(ptc.p())),
        ptc.frac());

    ptc_sec.set(species::electron);
    if (f) ptc_sec.set(f(ptc.template get<flagbits>(), species::electron));
    *(itr++) = ptc_sec;

    ptc_sec.set(species::positron);
    if (f) ptc_sec.set(f(ptc.template get<flagbits>(), species::positron));
    *(itr++) = std::move(ptc_sec);

  } else {
    *(itr++) = Particle<T, S>(ptc.q(), ptc.p() * (E_ph / apt::abs(ptc.p())),
                              ptc.frac(), species::photon,
                              f(ptc.template get<flagbits>(), species::photon));
  }
}

template <typename T, template <typename> class S>
void PhotonPairProduction(back_insert_iterator_t<T, S> itr, Ptc_t<T, S>& photon,
                          T, flagger_t f) {
  Particle<T, S> ptc_sec(
      photon.q(), photon.p() * std::sqrt(0.25 - 1.0 / apt::sqabs(photon.p())),
      photon.frac());
  ptc_sec.set(species::electron);
  if (f) ptc_sec.set(f(photon.template get<flagbits>(), species::electron));
  *(itr++) = ptc_sec;

  ptc_sec.set(species::positron);
  if (f) ptc_sec.set(f(photon.template get<flagbits>(), species::positron));
  *(itr++) = std::move(ptc_sec);

  photon.reset(flag::exist);  // void this photon
}
}  // namespace particle::scat

#include "particle/map.hpp"

namespace particle {
template <typename T, template <typename> class S>
map<Scat<T, S>> scat_map;

template <typename T, template <typename> class S>
void Scat<T, S>::Register(species sp) const {
  scat_map<T, S>.insert(sp, *this);
}

template <typename T, template <typename> class S>
void Scat<T, S>::Unregister(species sp) {
  scat_map<T, S>.erase(sp);
}

template <typename T, template <typename> class S>
Scat<T, S> Scat<T, S>::Get(species sp) {
  if (scat_map<T, S>.has(sp))
    return scat_map<T, S>[sp];
  else
    return {};
}
}  // namespace particle
