#ifndef _METRIC_CARTESIAN_HPP_
#define _METRIC_CARTESIAN_HPP_

namespace metric {
template <typename T>
struct Cartesian {
  template <int N>
  static inline T h(T = 0.0, T = 0.0, T = 0.0) noexcept {
    return 1.0;
  }

  template <int N>
  static inline T hh(T = 0.0, T = 0.0, T = 0.0) noexcept {
    static_assert(N == 0 || N == 1 || N == 2, "h index out of bounds");
    return 1.0;
  }

  static inline T hhh(T = 0.0, T = 0.0, T = 0.0) noexcept { return 1.0; }

  template <class X, class P>
  static inline auto geodesic_move(X& x, P& p, T dt, bool is_massive) noexcept {
    apt::array<T, P::NDim> dx;
    dt /= std::sqrt(is_massive + apt::sqabs(p));
    for (int i = 0; i < P::NDim; ++i) {
      dx[i] = p[i] * dt;
      x[i] += dx[i];
    }
    return dx;
  }
};
}  // namespace metric

#include "field/field.hpp"
#include "field/yee.hpp"

namespace field {
template <typename T>
constexpr T diff_one(T, T, T) noexcept {
  return 1.0;
}

// NOTE all derivatives are second order
template <int DGrid, typename T, int I, offset_t f_ofs>
constexpr T diff(const T& f, T, T, const apt::Grid<T, DGrid>& g,
                 const apt::Index<DGrid + 1>& s) noexcept {
  static_assert(DGrid > 1);
  static_assert(I >= 0 && I < 3);
  if constexpr (I >= DGrid) return 0.0;

  if constexpr (INSITU == f_ofs)
    return (*(&f + s[I]) - f) / g[I].delta();
  else
    return (f - *(&f - s[I])) / g[I].delta();
}

template <int DGrid, typename T, offset_t Ftype>
constexpr T (*Diff(int Fcomp,
                   int Icoord))(const T& f, T, T, const apt::Grid<T, DGrid>& g,
                                const apt::Index<DGrid + 1>& s) noexcept {
  if (0 == Fcomp) {
    if (1 == Icoord) return diff<DGrid, T, 1, !Ftype>;
    if (2 == Icoord) return diff<DGrid, T, 2, !Ftype>;
  } else if (1 == Fcomp) {
    if (2 == Icoord) return diff<DGrid, T, 2, !Ftype>;
    if (0 == Icoord) return diff<DGrid, T, 0, !Ftype>;
  } else {
    if (0 == Icoord) return diff<DGrid, T, 0, !Ftype>;
    if (1 == Icoord) return diff<DGrid, T, 1, !Ftype>;
  }
  return nullptr;
}
}  // namespace field

#endif
