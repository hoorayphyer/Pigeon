#ifndef  _APT_INTEGER_CLASS_HPP_
#define  _APT_INTEGER_CLASS_HPP_

#include <type_traits>

namespace apt {
  // TODO a safer way is to use compile type string const char* instead of int for LABEL, but this needs quite some trick
  template < typename T, int LABEL = 0 >
  struct Integer {
  public:
    using int_type = T;
    using self_type = Integer;
    constexpr Integer() noexcept = default;
    constexpr Integer(const Integer&) noexcept = default;
    constexpr Integer(Integer&&) noexcept = default;

    constexpr Integer( int_type a) noexcept : _n(a) {}

    constexpr operator int_type() const noexcept { return _n; }

    constexpr Integer& operator= ( const Integer& a ) noexcept = default;
    constexpr Integer& operator= ( Integer&& a ) noexcept = default;

    constexpr Integer& operator= ( int_type a ) noexcept {
      _n = a;
      return *this;
    }

    constexpr Integer& operator+=( int_type a ) noexcept {
      _n += a;
      return *this;
    }

    template < typename U = T >
    constexpr std::enable_if_t< std::is_signed_v<U>, Integer >
    operator-() const noexcept {
      return {-_n};
    }

  private:
    int_type _n{};
  };
}

namespace std {
  // Mimic enum behavior
  template < typename T, int L >
  struct underlying_type<apt::Integer<T,L>> {
    using type = typename apt::Integer<T,L>::int_type;
  };
}

#endif
