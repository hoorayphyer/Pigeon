#ifndef  _APT_VEC_HPP_
#define  _APT_VEC_HPP_

#include "apt/vec_expression.hpp"
#include "apt/array.hpp"

namespace apt {
  template < typename T, int N >
  struct Vec : public VecExpression<Vec<T,N>, T>,
               public VecModAssign< VecExpression<Vec<T,N>, T> > {
  private:
    array<T,N> _v{}; // NOTE {} here performs zero initialization.

  public:
    static constexpr auto NDim = N;
    using element_type = T;

    constexpr Vec() noexcept = default;

    template < typename... U, class = std::enable_if_t<sizeof...(U) == N> >
    constexpr Vec( U... args ) noexcept : _v{{ args... }} {}

    template < typename U >
    constexpr Vec( const array<U, N>& arr ) noexcept : _v{arr} {}

    template < typename E >
    constexpr Vec( const VecExpression<E>& vec ) noexcept {
      foreach<0,N>( [](auto& a, const auto& b){ a = b;}, _v, vec);
    }

    constexpr T operator[] ( int i ) const noexcept {
      return _v[i];
    }

    constexpr T& operator[] ( int i ) noexcept {
      return _v[i];
    }

  };
}

#endif
