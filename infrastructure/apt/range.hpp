#pragma once

#include "apt/array.hpp"
#include "apt/pair.hpp"

namespace apt {
struct Range {
 private:
  int _begin;
  int _end;
  apt::pair<int> _margin;

 public:
  Range() = default;
  Range(const Range&) = default;
  Range(Range&&) noexcept = default;
  Range& operator=(const Range&) = default;
  Range& operator=(Range&&) = default;

  constexpr Range(int begin, int end, apt::pair<int> margin) noexcept
      : _begin(begin), _end(end), _margin(margin) {}

  constexpr Range(int begin, int end) noexcept : Range(begin, end, {}) {}

  constexpr Range(int begin, int end, int margin) noexcept
      : Range(begin, end, {margin, margin}) {}

  constexpr bool operator==(const Range& other) const noexcept {
    return _begin == other._begin and _end == other._end and
           _margin[LFT] == other._margin[LFT] and
           _margin[RGT] == other._margin[RGT];
  }

  constexpr bool operator!=(const Range& other) const noexcept {
    return !(*this == other);
  }

  constexpr const int& begin() const noexcept { return _begin; }
  constexpr int& begin() noexcept { return _begin; }
  constexpr const int& end() const noexcept { return _end; }
  constexpr int& end() noexcept { return _end; }
  constexpr const auto& margin() const noexcept { return _margin; }
  constexpr auto& margin() noexcept { return _margin; }
  constexpr const auto& margin(bool i) const noexcept { return _margin[i]; }
  constexpr auto& margin(bool i) noexcept { return _margin[i]; }

  constexpr int size() const noexcept { return _end - _begin; }
  constexpr int far_begin() const noexcept { return _begin - _margin[LFT]; }
  constexpr int far_end() const noexcept { return _end + _margin[RGT]; }
  constexpr int full_size() const noexcept { return far_end() - far_begin(); }
};

template <int D>
constexpr array<Range, D> make_range(const array<int, D>& b,
                                     const array<int, D>& e,
                                     const array<pair<int>, D>& mg) noexcept {
  array<Range, D> res;
  for (int i = 0; i < D; ++i) res[i] = {b[i], e[i], mg[i]};
  return res;
}

template <int D>
constexpr array<Range, D> make_range(const array<int, D>& b,
                                     const array<int, D>& e, int mg) noexcept {
  array<Range, D> res;
  for (int i = 0; i < D; ++i) res[i] = {b[i], e[i], {mg, mg}};
  return res;
}
}  // namespace apt

#include "apt/array.hpp"
namespace apt::range {
template <int D>
constexpr array<int, D> begin(const array<Range, D>& range) noexcept {
  array<int, D> res;
  for (int i = 0; i < D; ++i) res[i] = range[i].begin();
  return res;
}

template <int D>
constexpr array<int, D> end(const array<Range, D>& range) noexcept {
  array<int, D> res;
  for (int i = 0; i < D; ++i) res[i] = range[i].end();
  return res;
}

template <int D>
constexpr array<pair<int>, D> margin(const array<Range, D>& range) noexcept {
  array<pair<int>, D> res;
  for (int i = 0; i < D; ++i) res[i] = range[i].margin();
  return res;
}

template <int D>
constexpr array<int, D> size(const array<Range, D>& range) noexcept {
  array<int, D> res;
  for (int i = 0; i < D; ++i) res[i] = range[i].size();
  return res;
}

template <int D>
constexpr array<int, D> far_begin(const array<Range, D>& range) noexcept {
  array<int, D> res;
  for (int i = 0; i < D; ++i) res[i] = range[i].far_begin();
  return res;
}

template <int D>
constexpr array<int, D> far_end(const array<Range, D>& range) noexcept {
  array<int, D> res;
  for (int i = 0; i < D; ++i) res[i] = range[i].far_end();
  return res;
}

template <int D>
constexpr array<int, D> full_size(const array<Range, D>& range) noexcept {
  array<int, D> res;
  for (int i = 0; i < D; ++i) res[i] = range[i].full_size();
  return res;
}

template <int D>
constexpr bool is_margin_uniform(const array<Range, D>& range) noexcept {
  int mg = range[0].margin()[LFT];
  for (int i = 0; i < D; ++i) {
    if (range[i].margin(LFT) != mg or range[i].margin(RGT) != mg) return false;
  }
  return true;
}

// setters
template <int D>
constexpr void set_begin(array<Range, D>& range,
                         const array<int, D>& beg) noexcept {
  for (int i = 0; i < D; ++i) range[i].begin() = beg[i];
}

template <int D>
constexpr void set_end(array<Range, D>& range,
                       const array<int, D>& end) noexcept {
  for (int i = 0; i < D; ++i) range[i].end() = end[i];
}

template <int D>
constexpr void set_margin(array<Range, D>& range,
                          const array<pair<int>, D>& margin) noexcept {
  for (int i = 0; i < D; ++i) range[i].margin() = margin[i];
}

template <int D>
constexpr void set_margin(array<Range, D>& range, int mg) noexcept {
  for (int i = 0; i < D; ++i) {
    range[i].margin(LFT) = mg;
    range[i].margin(RGT) = mg;
  }
}

// the following are handy when there is templated class deriving from it. See
// ActionBase or Haugbolle for example
template <int D>
constexpr const int& begin(const array<Range, D>& range, int i) noexcept {
  return range[i].begin();
}
template <int D>
constexpr int& begin(array<Range, D>& range, int i) noexcept {
  return range[i].begin();
}

template <int D>
constexpr const int& end(const array<Range, D>& range, int i) noexcept {
  return range[i].end();
}
template <int D>
constexpr int& end(array<Range, D>& range, int i) noexcept {
  return range[i].end();
}

template <int D>
constexpr const pair<int>& margin(const array<Range, D>& range,
                                  int i) noexcept {
  return range[i].margin();
}
template <int D>
constexpr pair<int>& margin(array<Range, D>& range, int i) noexcept {
  return range[i].margin();
}

template <int D>
constexpr int size(const array<Range, D>& range, int i) noexcept {
  return range[i].size();
}

template <int D>
constexpr int far_begin(const array<Range, D>& range, int i) noexcept {
  return range[i].far_begin();
}

template <int D>
constexpr int far_end(const array<Range, D>& range, int i) noexcept {
  return range[i].far_end();
}

template <int D>
constexpr int full_size(const array<Range, D>& range, int i) noexcept {
  return range[i].full_size();
}

// test empty. A range is empty when it's bulk is empty
constexpr bool is_empty(const Range& range) noexcept {
  return range.size() <= 0;
}
template <int D>
constexpr bool is_empty(const array<Range, D>& range) noexcept {
  for (int i = 0; i < D; ++i)
    if (is_empty(range[i])) return true;
  return false;
}
}  // namespace apt::range
