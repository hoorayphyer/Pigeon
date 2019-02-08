#ifndef  _APT_NUMERIC_HPP_
#define  _APT_NUMERIC_HPP_

#include "apt/vec_expression.hpp"

namespace apt {
  template <typename E1, typename E2>
  struct VecPlus : public VecExpression<VecPlus<E1,E2>> {
    static constexpr auto size = E1::size;
    const E1& e1;
    const E2& e2;

    template < int I >
    constexpr auto v() const noexcept {
      return e1.template v<I>() + e2.template v<I>();
    }
  };

  template <typename E1, typename E2>
  constexpr auto operator+(const E1& e1, const E2& e2) noexcept {
    return VecPlus<E1, E2>{e1, e2};
  }
}

namespace apt {
  template <typename E1, typename E2>
  struct VecMinus : public VecExpression<VecMinus<E1,E2>>{
    static constexpr auto size = E1::size;
    const E1& e1;
    const E2& e2;

    template < int I >
    constexpr auto v() const noexcept {
      return e1.template v<I>() - e2.template v<I>();
    }
  };

  template <typename E1, typename E2>
  constexpr auto operator-(const E1& e1, const E2& e2) noexcept {
    return VecMinus<E1, E2>{e1, e2};
  }
}

namespace apt {
  template <typename E1, typename T >
  struct VecTimes : public VecExpression<VecTimes<E1,T>> {
    static constexpr auto size = E1::size;
    const E1& e1;
    const T& t;

    template < int I >
    constexpr auto v() const noexcept {
      return e1.template v<I>() * t;
    }
  };

  template <typename E1, typename T>
  constexpr auto operator*(const E1& e1, const T& t) noexcept {
    return VecTimes<E1, T>{e1, t};
  }

  template <typename E1, typename T>
  constexpr auto operator/( const E1& e1, const T& t) noexcept {
    return VecTimes<E1, T>{e1, 1.0 / t};
  }
}

namespace apt {
  template < typename Vec1, typename Vec2 >
  constexpr lvec_ref<Vec1> operator += ( Vec1& v1, const Vec2& v2 ) noexcept {
    foreach<0,Vec1::size>([]( auto& a, const auto& b ) { a += b; }, v1, v2);
    return v1;
  }

  template < typename Vec1, typename U >
  constexpr lvec_ref<Vec1> operator *= ( Vec1& v1, const U& u ) noexcept {
    foreach<0,Vec1::size>([&u]( auto& a ) { a *= u; }, v1 );
    return v1;
  }
}

namespace apt {
  template <typename E1, typename E2>
  struct VecCross : public VecExpression<VecCross<E1,E2>> {
    static_assert( E1::size == 3 && E2::size == 3 );
    static constexpr auto size = E1::size;

    const E1& e1;
    const E2& e2;

    template < int I >
    constexpr auto v() const noexcept {
      return e1.template v<(I+1)%3>() * e2.template v<(I+2)%3>()
        - e1.template v<(I+2)%3>() * e2.template v<(I+1)%3>();
    }
  };

  template < typename E1, typename E2 >
  constexpr auto cross( const E1& e1, const E2& e2 ) noexcept {
    return VecCross<E1,E2>{e1, e2};
  }
}

#include <cmath>
namespace apt {
  template < int I, typename E1, typename E2 >
  constexpr auto dot ( const VecExpression<E1>& v1, const VecExpression<E2>& v2 ) noexcept {
    if constexpr ( I == E1::size ) return 0;
    else return v1.template v<I>() * v2.template v<I>() + dot<I+1>(v1, v2);
  }

  template < typename E1, typename E2 >
  constexpr auto dot ( const VecExpression<E1>& v1, const VecExpression<E2>& v2 ) noexcept {
    return dot<0,E1,E2>( v1, v2 );
  }

  template < typename E >
  constexpr auto sqabs( const VecExpression<E>& vec ) noexcept {
    return dot(vec, vec);
  }

  template < typename E >
  constexpr auto abs( const VecExpression<E>& vec ) noexcept {
    if constexpr ( E::size == 1 ) return std::abs( vec.template v<0>() );
    else return std::sqrt( sqabs(vec) );
  }

}


#endif
