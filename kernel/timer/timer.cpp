#include "timer/timer.hpp"
#include <chrono>
using namespace std::chrono;

namespace tmr {
  inline auto now() {
    return high_resolution_clock::now();
  }
  using timepoint_t = decltype(now());

  TimeStamp::TimeStamp()
    : _t( now() ) {}

  double TimeStamp::lapse() {
    auto dur = duration_cast<nanoseconds>( now() - std::any_cast<timepoint_t>(_t) );
    return static_cast<double>( dur.count() ) / 1e6;
  }

  double TimeStamp::operator- ( const TimeStamp& other ) {
    auto dur = duration_cast<nanoseconds>( std::any_cast<timepoint_t>(_t) - std::any_cast<timepoint_t>(other._t) );
    return static_cast<double>( dur.count() ) / 1e6;
  }

}
