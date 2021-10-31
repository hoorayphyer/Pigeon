#pragma once

#include "particle/state.hpp"

namespace particle {
template <typename Ptc, typename Vec = typename Ptc::vec_type,
          typename T = typename Vec::element_type,
          typename State = typename Ptc::state_type>
struct PtcExpression : public particle::StateExpression<Ptc, State> {
  static constexpr int NDim = Ptc::NDim;
  using vec_type = Vec;
  using state_type = State;

  constexpr Vec& q() noexcept { return static_cast<Ptc&>(*this).q(); }
  constexpr const Vec& q() const noexcept {
    return static_cast<const Ptc&>(*this).q();
  }
  constexpr T& q(int i) noexcept { return static_cast<Ptc&>(*this).q(i); }
  constexpr const T& q(int i) const noexcept {
    return static_cast<const Ptc&>(*this).q(i);
  }

  constexpr Vec& p() noexcept { return static_cast<Ptc&>(*this).p(); }
  constexpr const Vec& p() const noexcept {
    return static_cast<const Ptc&>(*this).p();
  }
  constexpr T& p(int i) noexcept { return static_cast<Ptc&>(*this).p(i); }
  constexpr const T& p(int i) const noexcept {
    return static_cast<const Ptc&>(*this).p(i);
  }

  constexpr T& frac() noexcept { return static_cast<Ptc&>(*this).frac(); }
  constexpr const T& frac() const noexcept {
    return static_cast<const Ptc&>(*this).frac();
  }

  // NOTE state() are defined in StateExpression
};

}  // namespace particle
