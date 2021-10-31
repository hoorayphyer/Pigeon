#pragma once

#include "apt/array.hpp"

namespace apt {
constexpr int pow3(int i) noexcept {
  if (0 == i)
    return 1;
  else
    return 3 * pow3(i - 1);
}
}  // namespace apt

namespace apt {
constexpr int L = 0;
constexpr int C = 1;
constexpr int R = 2;
}  // namespace apt
