#ifndef  _APT_VEC_HPP_
#define  _APT_VEC_HPP_

#include "apt/vec_expression.hpp"
#include <array>

namespace apt {
  template < typename T, int N >
  struct Vec : public VecExpression<Vec<T,N>, T> {
  private:
    std::array<T,N> _v;

  public:
    static constexpr auto size = N;
    using value_type = T;

    constexpr Vec() = default;

    template < typename U >
    constexpr Vec( U x[N] ) noexcept {
      foreach<0,N>( [](auto& a, const auto& b){ a = b;}, _v, x);
    }

    template < typename E >
    constexpr Vec( const VecExpression<E>& vec ) noexcept {
      foreach<0,N>( [](auto& a, const auto& b){ a = b;}, _v, vec);
    }

    template < int I >
    constexpr T v() const noexcept {
      return std::get<I>(_v);
    }

    template < int I >
    constexpr T& v() noexcept {
      return std::get<I>(_v);
    }

  };
}

#endif
