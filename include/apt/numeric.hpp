#ifndef  _APT_NUMERIC_HPP_
#define  _APT_NUMERIC_HPP_

#include "apt/vec_expression.hpp"

namespace apt {
  template <typename E1, typename E2, typename T = typename E1::value_type>
  struct VecPlus : public VecExpression<VecPlus<E1,E2,T>, T> {
  private:
    const E1& _e1;
    const E2& _e2;
  public:
    using value_type = T;
    static constexpr int size = E1::size;

    constexpr VecPlus( const E1& e1, const E2& e2 ) noexcept
      : _e1(e1), _e2(e2){}

    template < int I >
    constexpr T v() const noexcept {
      return _e1.template v<I>() + _e2.template v<I>();
    }
  };

  template <typename E1, typename E2, typename T = typename E1::value_type>
  struct VecMinus : public VecExpression<VecMinus<E1,E2,T>, T>{
  private:
    const E1& _e1;
    const E2& _e2;

  public:
    using value_type = T;
    static constexpr int size = E1::size;

    constexpr VecMinus( const E1& e1, const E2& e2 ) noexcept
      : _e1(e1), _e2(e2){}

    template < int I >
    constexpr T v() const noexcept {
      return _e1.template v<I>() - _e2.template v<I>();
    }
  };

  template <typename E1, typename T, typename U = typename E1::value_type >
  struct VecTimes : public VecExpression<VecTimes<E1,T>, U> {
  private:
    const E1& _e1;
    const T& _t;

  public:
    using value_type = T;
    static constexpr int size = E1::size;

    constexpr VecTimes( const E1& e1, const T& t ) noexcept
      : _e1(e1), _t(t) {}

    template < int I >
    constexpr U v() const noexcept {
      return _e1.template v<I>() * _t;
    }
  };
}

namespace {
  template <typename E1, typename E2>
  constexpr auto operator+(const apt::VecExpression<E1>& e1,
                           const apt::VecExpression<E2>& e2) noexcept {
    return apt::VecPlus<apt::VecExpression<E1>, apt::VecExpression<E2>>(e1, e2);
  }

  template <typename E1, typename E2>
  constexpr auto operator-(const apt::VecExpression<E1>& e1,
                           const apt::VecExpression<E2>& e2) noexcept {
    return apt::VecMinus<apt::VecExpression<E1>, apt::VecExpression<E2>>(e1, e2);
  }

  template <typename E1, typename T>
  constexpr auto operator*(const apt::VecExpression<E1>& e1, const T& t) noexcept {
    return apt::VecTimes<apt::VecExpression<E1>, T>(e1, t);
  }

  template <typename E1, typename T>
  constexpr auto operator/( const apt::VecExpression<E1>& e1, const T& t) noexcept {
    return apt::VecTimes<E1, T>(e1, 1.0 / t);
  }

  template <typename E1, typename E2 >
  constexpr auto& operator+= ( apt::VecExpression<E1>& v1, const apt::VecExpression<E2>& v2 ) noexcept {
    apt::foreach<0, E1::size>
      ( []( auto& x, const auto& y ) noexcept { x += y; }, v1, v2 );
    return v1;
  }

  template <typename E1, typename E2 >
  constexpr auto& operator-= ( apt::VecExpression<E1>& v1, const apt::VecExpression<E2>& v2 ) noexcept {
    apt::foreach<0, E1::size>
      ( []( auto& x, const auto& y ) noexcept { x -= y; }, v1, v2 );
    return v1;
  }

  template <typename E, typename T = typename E::value_type>
  constexpr auto& operator*= ( apt::VecExpression<E>& v, const T& t ) noexcept {
    apt::foreach<0, E::size>
      ( [&t]( auto& x ) noexcept { x *= t; }, v );
    return v;
  }

  template <typename E, typename T = typename E::value_type>
  constexpr auto& operator/= ( apt::VecExpression<E>& v, const T& t ) noexcept {
    apt::foreach<0, E::size>
      ( [&t]( auto& x ) noexcept { x /= t; }, v );
    return v;
  }
}

namespace apt {
  template <typename E1, typename E2, typename T = typename E1::value_type>
  struct VecCross : public VecExpression<VecCross<E1,E2,T>, T> {
  private:
    const E1& _e1;
    const E2& _e2;

  public:
    using value_type = T;
    static constexpr int size = E1::size;

    constexpr VecCross( const E1& e1, const E2& e2 ) noexcept
      : _e1(e1), _e2(e2) {}

    template < int I >
    constexpr T v() const noexcept {
      return _e1.template v<(I+1)%3>() * _e2.template v<(I+2)%3>()
        - _e1.template v<(I+2)%3>() * _e2.template v<(I+1)%3>();
    }
  };

  template < typename E1, typename E2 >
  constexpr auto cross( const E1& e1, const E2& e2 ) noexcept {
    return VecCross<E1,E2>(e1, e2);
  }
}

#include <cmath>
namespace apt {
  template < int I, typename E1, typename E2 >
  constexpr auto dot ( const VecExpression<E1>& v1, const VecExpression<E2>& v2 ) noexcept {
  }

  template < typename E1, typename E2, int I = 0 >
  constexpr auto dot ( const VecExpression<E1>& v1, const VecExpression<E2>& v2 ) noexcept {
    constexpr int END = E1::size;
    if constexpr ( I == END ) return 0;
    else return v1.template v<I>() * v2.template v<I>() + dot<E1, E2, I+1>(v1, v2);
  }

  template < typename E >
  constexpr auto sqabs( const VecExpression<E>& vec ) noexcept {
    return dot(vec, vec);
  }

  template < typename E >
  constexpr auto abs( const VecExpression<E>& vec ) noexcept {
    constexpr int N = E::size;
    if constexpr ( N == 1 ) return std::abs( vec.template v<0>() );
    else return std::sqrt( sqabs(vec) );
  }

}


#endif
