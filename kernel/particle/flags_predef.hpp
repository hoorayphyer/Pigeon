#pragma once

namespace particle {
// 16 flags
enum class flag : unsigned int {
  exist = 0,
  secondary,
  ignore_force,
  _3,
  _4,
  _5,
  _6,
  _7,
  _8,
  _9,
  _10,
  _11,
  _12,
  _13,
  _14,
  _15
};

constexpr bool is_reserved_flag(flag f) noexcept {
  return static_cast<unsigned int>(f) < 3;
}
}  // namespace particle
