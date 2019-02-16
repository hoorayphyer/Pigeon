#ifndef _APT_VEC_SLICE_HPP_
#define _APT_VEC_SLICE_HPP_

#include "apt/vec_expression.hpp"

namespace apt {
  template < int Begin, int End, typename E, typename T = typename E::value_type >
  struct VecSlice : public VecExpression<VecSlice<Begin, End, E, T>, T> {
  private:
    static_assert(0 <= Begin && Begin <= End);
    E& _e;

  public:
    static constexpr int size = End - Begin;
    using value_type = T;

    constexpr VecSlice( E& e ) noexcept : _e(e) {}

    template < int I >
    constexpr std::enable_if_t< (size > 0), T >
    v() const noexcept {
      return _e.template v< Begin + I >();
    }

    template < int I >
    constexpr std::enable_if_t< (size > 0), T >&
    v() noexcept {
      return _e.template v< Begin + I >();
    }

  };
}


#endif
