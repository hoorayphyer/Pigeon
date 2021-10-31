#pragma once
#include "apt/foreach.hpp"

namespace apt {
template <typename E, typename Real = typename E::element_type>
class VecExpression {
 public:
  static constexpr int NDim = E::NDim;
  using element_type = Real;

  constexpr Real operator[](int i) const noexcept {
    return static_cast<const E&>(*this)[i];
  }
};
}  // namespace apt

namespace apt {
template <typename Vec>
struct VecModAssign;

template <typename E1, typename Real>
struct VecModAssign<VecExpression<E1, Real> > {
 private:
  constexpr auto& self() noexcept { return static_cast<E1&>(*this); };

 public:
  template <typename E2>
  constexpr E1& operator+=(const apt::VecExpression<E2>& v2) noexcept {
    apt::foreach<0, E1::NDim>([](auto& x, const auto& y) noexcept { x += y; },
                              self(), v2);
    return self();
  }

  template <typename E2>
  constexpr E1& operator-=(const apt::VecExpression<E2>& v2) noexcept {
    apt::foreach<0, E1::NDim>([](auto& x, const auto& y) noexcept { x -= y; },
                              self(), v2);
    return self();
  }

  constexpr E1& operator*=(Real t) noexcept {
    apt::foreach<0, E1::NDim>([&t](auto& x) noexcept { x *= t; }, self());
    return self();
  }

  constexpr E1& operator/=(Real t) noexcept {
    apt::foreach<0, E1::NDim>([&t](auto& x) noexcept { x /= t; }, self());
    return self();
  }
};

}  // namespace apt
