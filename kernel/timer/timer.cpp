#include "timer/timer.hpp"
#include <chrono>
using namespace std::chrono;

namespace tmr {
  inline auto now() {
    return high_resolution_clock::now();
  }

  template < typename T >
  constexpr auto convert( const T& reading_in_ns, std::string unit, std::string unit_fallback ) noexcept {
    if ( "ms" == unit ) return reading_in_ns / 1e6;
    if ( "us" == unit ) return reading_in_ns / 1e3;
    if ( "ns" == unit ) return reading_in_ns;
    if ( "s" == unit ) return  reading_in_ns/ 1e9;
    return convert( reading_in_ns, unit_fallback, "ms"); // return default unit if unspecified unit is passed in
  }

  using timepoint_t = decltype(now());

  Timestamp::Timestamp()
    : _t( now() ) {}

  double Timestamp::lapse(std::string unit) const {
    auto dur = duration_cast<nanoseconds>( now() - std::any_cast<timepoint_t>(_t) );
    return convert( static_cast<double>( dur.count() ), unit, _unit );
  }

  double Timestamp::operator- ( const Timestamp& other ) {
    auto dur = duration_cast<nanoseconds>( std::any_cast<timepoint_t>(_t) - std::any_cast<timepoint_t>(other._t) );
    return convert( static_cast<double>( dur.count() ), _unit, "ms" );
  }

}
