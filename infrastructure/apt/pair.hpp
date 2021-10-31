#pragma once

namespace apt {
template <typename T>
struct pair {
  T lft;
  T rgt;

  template <typename U>
  constexpr pair& operator=(const apt::pair<U>& b) noexcept {
    lft = b.lft;
    rgt = b.rgt;
    return *this;
  }

  constexpr const T& operator[](bool lr) const noexcept {
    return lr ? rgt : lft;
  }
  constexpr T& operator[](bool lr) noexcept { return lr ? rgt : lft; }
};

template <typename T>
pair<T&> tie(T& lft, T& rgt) noexcept {
  return {lft, rgt};
}
}  // namespace apt

inline constexpr bool LFT = false;
inline constexpr bool RGT = true;
