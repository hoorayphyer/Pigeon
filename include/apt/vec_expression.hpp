#ifndef _APT_VEC_EXPRESSION_HPP_
#define _APT_VEC_EXPRESSION_HPP_

namespace apt {

  template <typename E,
            typename T = typename E::element_type,
            bool is_lvalue = E::is_lvalue>
  class VecExpression;

  template <typename E, typename T>
  class VecExpression<E,T,false> {
  public:
    static constexpr int NDim = E::NDim;
    using element_type = T;
    static constexpr bool is_lvalue = false;

    constexpr T operator[] ( int i ) const noexcept {
      return static_cast<const E&>(*this)[i];
    }

  };

  template <typename E, typename T>
  class VecExpression<E,T,true> {
  public:
    static constexpr int NDim = E::NDim;
    using element_type = T;
    static constexpr bool is_lvalue = true;

    constexpr T operator[] ( int i ) const noexcept {
      return static_cast<const E&>(*this)[i];
    }

    constexpr T& operator[] ( int i ) noexcept {
      return static_cast<E&>(*this)[i];
    }
  };
}

#include "apt/foreach.hpp"

namespace std{
  // TODO can this work on vVec?
  template < typename E1, typename T1, typename E2, typename T2 >
  void swap ( apt::VecExpression<E1,T1>& a, apt::VecExpression<E2,T2>& b ) noexcept {
    apt::foreach<0, E1::NDim>
      ( []( auto& x, auto& y ) noexcept { std::swap(x,y); }, a, b );
  }
}

#endif
