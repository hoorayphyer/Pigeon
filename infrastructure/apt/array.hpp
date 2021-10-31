#pragma once

namespace apt {
template <typename T, int D>
struct array {
  T _data[D]{};

  using element_type = T;
  static constexpr int NDim = D;

  // TODO bound checks.
  constexpr const T& operator[](int i) const noexcept { return _data[i]; }
  constexpr T& operator[](int i) noexcept { return _data[i]; }

  constexpr bool operator!=(const array& other) const noexcept {
    for (int i = 0; i < D; ++i) {
      if (_data[i] != other[i]) return true;
    }
    return false;
  }

  constexpr bool operator==(const array& other) const noexcept {
    return !(*this != other);
  }

  constexpr T* begin() noexcept { return _data; }
  constexpr const T* begin() const noexcept { return _data; }

  constexpr const T* end() const noexcept { return _data + D; }

  constexpr const T& back() const noexcept { return _data[D - 1]; }
  constexpr T& back() noexcept { return _data[D - 1]; }

  constexpr int size() const noexcept { return NDim; }
};

// a C-array is not allowed to have size 0
template <typename T>
struct array<T, 0> {
  static constexpr int NDim = 0;
};
}  // namespace apt
