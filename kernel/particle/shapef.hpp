#pragma once

#include <cmath>

// value here denotes length of support
namespace particle {
enum class shape : int {
  Nearest_Grid_Point = 1,
  Cloud_In_Cell = 2,
  Triangular_Cloud = 3,
  Piecewise_Cubic_Spline = 4
};
}

namespace particle {
template <shape S>
struct shapef_t {
  constexpr static int support() noexcept { return static_cast<int>(S); }

  template <typename T>
  constexpr T operator()(T dx) const noexcept {
    if constexpr (shape::Nearest_Grid_Point == S) {
      return static_cast<T>(dx < 0.5 && dx >= -0.5);
    } else {
      dx = std::abs(dx);
      if constexpr (shape::Cloud_In_Cell == S) {
        return std::max<T>(1.0 - dx, 0.0);
      } else if (shape::Triangular_Cloud == S) {
        return static_cast<T>(dx < 0.5) * (0.75 - dx * dx) +
               static_cast<T>(dx >= 0.5 && dx < 1.5) * 0.5 * (1.5 - dx) *
                   (1.5 - dx);
      } else if (shape::Piecewise_Cubic_Spline == S) {
        return static_cast<T>(dx < 1.0) *
                   (2.0 / 3.0 - dx * dx * (1.0 - 0.5 * dx)) +
               static_cast<T>(dx >= 1.0 && dx < 2.0) * (2.0 - dx) * (2.0 - dx) *
                   (2.0 - dx) / 6.0;
      }
    }
  }
};

template <class SF, int D>
struct induced_shapef_t {
  static_assert(D > 0);
  constexpr static int support() noexcept { return SF::support() / D; }

  // For fields defined at MIDWAY
  template <typename T>
  constexpr T operator()(T dx) const noexcept {
    if constexpr (D == 1)
      return SF()(dx);
    else {
      static_assert(D % 2 == 0);
      SF sf;
      T sum = 0;
      for (int i = 0; i < D; ++i) sum += sf(D * dx + D / 2 - i - 0.5);
      return sum;
    }
  }
};

}  // namespace particle
