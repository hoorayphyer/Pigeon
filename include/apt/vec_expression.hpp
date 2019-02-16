#ifndef _APT_VEC_EXPRESSION_HPP_
#define _APT_VEC_EXPRESSION_HPP_

namespace apt {
  template <typename E, typename T = typename E::value_type>
  class VecExpression {
  public:
    static constexpr int size = E::size;
    using value_type = T;

    // NOTE .template v<I>() is necessary to tell compiler that v is a template. It works like typename
    template < int I >
    constexpr T v() const noexcept {
      return static_cast<const E&>(*this).template v<I>();
    }

    template < int I >
    constexpr T& v() noexcept { // const because it won't change what references are stored
      return static_cast<E&>(*this).template v<I>();
    }
  };
}

namespace std {
  // define this so as to be used in apt::foreach
  template < int I, typename E >
  constexpr auto get( const apt::VecExpression<E>& vec ) noexcept {
    return vec.template v<I>();
  }

  template < int I, typename E >
  constexpr auto& get( apt::VecExpression<E>& vec ) noexcept {
    return vec.template v<I>();
  }
}

#include "apt/foreach.hpp"

namespace std{
  // TODO can this work on vVec?
  template < typename E1, typename T1, typename E2, typename T2 >
  void swap ( apt::VecExpression<E1,T1>& a, apt::VecExpression<E2,T2>& b ) noexcept {
    apt::foreach<0, E1::size>
      ( []( auto& x, auto& y ) noexcept { std::swap(x,y); }, a, b );
  }
}

#endif
