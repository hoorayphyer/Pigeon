#pragma once

// TODOL all the stuff under this {} are meant to be user-specified. Here the
// pulsar in LogSpherical is used
namespace particle::force {
// LogSpherical
template <typename T, template <typename> class Specs,
          template <typename, template <typename> class> class Ptc_t>
void landau0(Ptc_t<T, Specs>& ptc, T dt, const apt::Vec<T, Specs<T>::Dim>& E,
             const apt::Vec<T, Specs<T>::Dim>& B, T B_thr) noexcept {
  using Vec = apt::Vec<T, Specs<T>::Dim>;
  if (apt::sqabs(B) < B_thr * B_thr) return;

  auto EB2 = apt::dot(E, B);
  EB2 = EB2 * EB2;
  auto B2_E2 = apt::sqabs(B) - apt::sqabs(E);
  // calculate E'^2
  auto Ep2 = 2 * EB2 / (std::sqrt(B2_E2 * B2_E2 + 4 * EB2) + B2_E2);
  Vec beta_ExB = apt::cross(E, B) / (apt::sqabs(B) + Ep2);
  // find B' modulo gamma_ExB
  Vec Bp = B - apt::cross(beta_ExB, E);
  // obtain the momentum with perpendicular components damped
  ptc.p() = Bp * (apt::dot(ptc.p(), Bp) / apt::sqabs(Bp));
  ptc.p() += beta_ExB * std::sqrt((1.0 + apt::sqabs(ptc.p())) /
                                  (1.0 - apt::sqabs(beta_ExB)));
}
}  // namespace particle::force
