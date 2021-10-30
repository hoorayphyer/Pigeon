#include "timer.hpp"

#include <chrono>
using namespace std::chrono;

namespace tmr {
inline auto now() { return high_resolution_clock::now(); }

TDur Duration::val(std::string unit) const noexcept {
  if ("ms" == unit) return _val / 1e6;
  if ("us" == unit) return _val / 1e3;
  if ("ns" == unit) return _val;
  if ("s" == unit) return _val / 1e9;
  return val(_unit);  // return default unit if unspecified unit is passed in
}

using timepoint_t = decltype(now());

Timestamp::Timestamp() : _t(now()) {}

Duration Timestamp::lapse() const {
  auto dur = duration_cast<nanoseconds>(now() - std::any_cast<timepoint_t>(_t));
  return {static_cast<TDur>(dur.count())};
}

Duration Timestamp::operator-(const Timestamp& other) {
  auto dur = duration_cast<nanoseconds>(std::any_cast<timepoint_t>(_t) -
                                        std::any_cast<timepoint_t>(other._t));
  return {static_cast<TDur>(dur.count())};
}

}  // namespace tmr
