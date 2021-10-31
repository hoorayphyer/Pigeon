#pragma once

#include <cstddef>

namespace apt {
// convention: Pos counts from the right. (Pos+N, Pos].
// NOTE bit shift operator << returns of type of the left operand. So if LHS has
// fewer bits than RHS, the return result will simply be zero, i.e., information
// is lost, hence all the casts
template <std::size_t Pos, std::size_t N, typename T>
constexpr T getbits(const T &x) noexcept {
  return static_cast<T>((x >> Pos) & ~(~static_cast<std::size_t>(0) << N));
}

template <std::size_t Pos, std::size_t N, typename U, typename T>
constexpr void setbits(T &x, U y) noexcept {
  // NOTE assume y has all nonzero bits in the last N bits
  x &= ~(~(~static_cast<std::size_t>(0) << N) << Pos);
  x |= ((static_cast<std::size_t>(y) & ~(~static_cast<std::size_t>(0) << N))
        << Pos);
}
}  // namespace apt
