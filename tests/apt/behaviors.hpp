#include <experimental/type_traits>

template <template < class > class Behavior, typename Vec>
constexpr bool is( const Vec& ) noexcept {
  return std::experimental::is_detected_v<Behavior,Vec>;
}

namespace bhv {
  template <typename Vec>
  using lvec = decltype( std::get<0>(std::declval<Vec>()) = 0.0 );

  // template <typename Vec>
  // using vt_able = decltype( std::tie( std::get<0>(std::declval<Vec>()) ) );
}
