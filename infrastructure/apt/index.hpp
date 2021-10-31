#pragma once

#include "apt/array.hpp"

namespace apt {
template <int D>
using Index = array<int, D>;
}

template <int D>
constexpr apt::Index<D> operator+(const apt::Index<D>& ind_a,
                                  const apt::Index<D>& ind_b) noexcept {
  apt::Index<D> res;
  for (int i = 0; i < D; ++i) res[i] = ind_a[i] + ind_b[i];
  return res;
}

template <int D>
constexpr apt::Index<D> operator-(const apt::Index<D>& ind_a,
                                  const apt::Index<D>& ind_b) noexcept {
  apt::Index<D> res;
  for (int i = 0; i < D; ++i) res[i] = ind_a[i] - ind_b[i];
  return res;
}

template <int D>
constexpr apt::Index<D>& operator+=(apt::Index<D>& ind_a,
                                    const apt::Index<D>& ind_b) noexcept {
  for (int i = 0; i < D; ++i) ind_a[i] += ind_b[i];
  return ind_a;
}

template <int D>
constexpr bool operator==(const apt::Index<D>& ind_a,
                          const apt::Index<D>& ind_b) noexcept {
  bool res = true;
  for (int i = 0; i < D; ++i) res = (res and (ind_a[i] == ind_b[i]));
  return res;
}

template <int D>
constexpr bool operator!=(const apt::Index<D>& ind_a,
                          const apt::Index<D>& ind_b) noexcept {
  return !(ind_a == ind_b);
}

namespace apt {
struct Longidx {  // longitudianl index
 private:
  const int _dir = 0;
  int _val = 0;

 public:
  Longidx() = default;
  Longidx(const Longidx&) = default;
  Longidx(Longidx&&) noexcept = default;

  constexpr Longidx& operator=(const Longidx&) = default;
  constexpr Longidx& operator=(Longidx&&) noexcept = default;
  constexpr Longidx& operator=(int a) noexcept {
    _val = a;
    return *this;
  }

  constexpr Longidx(int dir, int val) noexcept : _dir(dir), _val(val) {}

  constexpr int dir() const noexcept { return _dir; }
  constexpr int val() const noexcept { return _val; }

  template <int D>
  constexpr Index<D> operator+(Index<D> x) const noexcept {
    x[_dir] += _val;
    return x;
  }

  template <int D>
  constexpr Index<D> operator-(Index<D> x) const noexcept {
    for (int i = 0; i < D; ++i) x[i] *= -1;
    x[_dir] += _val;
    return x;
  }

  constexpr Longidx operator+(int a) const noexcept {
    auto res(*this);
    res._val += a;
    return res;
  }
  constexpr Longidx operator-(int a) const noexcept {
    auto res(*this);
    res._val -= a;
    return res;
  }
  constexpr Longidx operator*(int a) const noexcept {
    auto res(*this);
    res._val *= a;
    return res;
  }
  constexpr Longidx operator/(int a) const noexcept {
    auto res(*this);
    res._val /= a;
    return res;
  }

  constexpr Longidx& operator+=(int a) noexcept {
    _val += a;
    return *this;
  }
  constexpr Longidx& operator-=(int a) noexcept {
    _val -= a;
    return *this;
  }
  constexpr Longidx& operator*=(int a) noexcept {
    _val *= a;
    return *this;
  }
  constexpr Longidx& operator/=(int a) noexcept {
    _val /= a;
    return *this;
  }

  constexpr Longidx& operator++() noexcept {
    ++_val;
    return *this;
  }
  constexpr Longidx& operator--() noexcept {
    --_val;
    return *this;
  }
  constexpr Longidx operator++(int) noexcept {
    auto res = *this;
    ++_val;
    return res;
  }
  constexpr Longidx operator--(int) noexcept {
    auto res = *this;
    --_val;
    return res;
  }

  constexpr bool operator<(int a) const noexcept { return _val < a; }
  constexpr bool operator<=(int a) const noexcept { return _val <= a; }
  constexpr bool operator>(int a) const noexcept { return _val > a; }
  constexpr bool operator>=(int a) const noexcept { return _val >= a; }
  constexpr bool operator!=(int a) const noexcept { return _val != a; }
  constexpr bool operator==(int a) const noexcept { return _val == a; }
};
}  // namespace apt

template <int D>
constexpr apt::Index<D> operator+(apt::Index<D> ind,
                                  const apt::Longidx& l) noexcept {
  ind[l.dir()] += l.val();
  return ind;
}

template <int D>
constexpr apt::Index<D> operator-(apt::Index<D> ind,
                                  const apt::Longidx& l) noexcept {
  ind[l.dir()] -= l.val();
  return ind;
}

constexpr apt::Longidx operator+(int a, apt::Longidx l) noexcept {
  l += a;
  return l;
}

constexpr apt::Longidx operator*(int a, apt::Longidx l) noexcept {
  l *= a;
  return l;
}

constexpr apt::Longidx operator-(int a, apt::Longidx l) noexcept {
  l = a - l.val();
  return l;
}
