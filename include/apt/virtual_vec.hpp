#ifndef  _APT_VIRTUAL_VEC_HPP_
#define  _APT_VIRTUAL_VEC_HPP_

#include "apt/vec_expression.hpp"
#include <tuple>

// NOTE virtual vector, or a vector proxy
namespace apt {
  namespace impl {
    template < typename T, int N > struct ref_tuple;
    template < typename T > struct ref_tuple<T,0> { using type = std::tuple<>; };
    template < typename T > struct ref_tuple<T,1> { using type = std::tuple<T&>; };
    template < typename T > struct ref_tuple<T,2> { using type = std::tuple<T&,T&>; };
    template < typename T > struct ref_tuple<T,3> { using type = std::tuple<T&,T&,T&>; };
    template < typename T, int N >
    using ref_tuple_t = typename ref_tuple<T,N>::type;
  }

  template < typename T, int N >
  struct vVec : public VecExpression<vVec<T,N>, T, true> {
  private:
    using tuple_type = impl::ref_tuple_t<T,N>;
    tuple_type _v;

    template < typename E, std::size_t... I > // NOTE using size_t instead of int is to make compiling work
    constexpr vVec( VecExpression<E>&& vec, std::index_sequence<I...> ) noexcept
      : _v( std::get<I>(vec)... ) {}

    template < std::size_t... I >
    constexpr vVec( std::array<T,N>& arr, std::index_sequence<I...> ) noexcept
      : vVec(std::get<I>(arr)...) {}

  public:
    using value_type = T;
    static constexpr int size = N;

    template < int I >
    constexpr T v() const noexcept {
      return std::get<I>(_v);
    }

    template < int I >
    constexpr T& v() noexcept {
      return std::get<I>(_v);
    }

    template < typename... U >
    constexpr vVec( U&... u ) noexcept
      : _v(u...) {};

    constexpr vVec( std::array<T,N>& arr ) noexcept
      : vVec( arr, std::make_index_sequence<N>{} ) {}

    constexpr vVec( tuple_type&& tpl ) noexcept
      : _v(std::move(tpl)) {};

    constexpr vVec( vVec&& vec ) noexcept = default;

    template <typename E>
    constexpr vVec( VecExpression<E>&& vec ) noexcept
      : vVec( std::move(vec), std::make_index_sequence<E::size>{} ) {}

    vVec() = delete;
    vVec( const vVec& ) = delete; // because it breaks copy sematics
    template <typename E>
    vVec( const VecExpression<E>& ) = delete;
  };

}

#endif
