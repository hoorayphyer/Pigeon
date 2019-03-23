#ifndef  _APT_VIRTUAL_VEC_HPP_
#define  _APT_VIRTUAL_VEC_HPP_

#include "apt/vec_expression.hpp"
#include "apt/array.hpp"
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

    template < typename E, std::size_t... I >
    constexpr vVec( VecExpression<E>&& vec, std::index_sequence<I...> ) noexcept
      : _v( vec[I]... ) {}

    template < std::size_t... I >
    constexpr vVec( array<T,N>& arr, std::index_sequence<I...> ) noexcept
      : vVec(arr[I]...) {}

  public:
    using element_type = T;
    static constexpr int NDim = N;

    constexpr T operator[] ( int i ) const noexcept {
      if ( 0 == i ) return std::get<0>(_v);
      if constexpr ( NDim > 1 )
                     if ( 1 == i ) return std::get<1>(_v);
      if constexpr ( NDim > 2 )
                     if ( 2 == i ) return std::get<2>(_v);
    }

    constexpr T& operator[] ( int i ) noexcept {
      if ( 0 == i ) return std::get<0>(_v);
      if constexpr ( NDim > 1 )
                     if ( 1 == i ) return std::get<1>(_v);
      if constexpr ( NDim > 2 )
                     if ( 2 == i ) return std::get<2>(_v);
    }

    template < typename... U >
    constexpr vVec( U&... u ) noexcept
      : _v(std::tie(u...)) {};

    constexpr vVec( array<T,N>& arr ) noexcept
      : vVec( arr, std::make_index_sequence<N>{} ) {}

    constexpr vVec( tuple_type&& tpl ) noexcept
      : _v(std::move(tpl)) {};

    constexpr vVec( vVec&& vec ) noexcept = default;

      // : vVec( std::move(vec), std::make_index_sequence<E::NDim>{} ) {}

    vVec() = delete;
    vVec( const vVec& ) = delete; // because it breaks copy sematics

    constexpr vVec& operator= ( const vVec& vec ) noexcept {
      _v = vec._v;
      return *this;
    }

    template <typename E>
    constexpr vVec& operator= ( const VecExpression<E>& vec ) noexcept {
      std::get<0>(_v) = vec[0];
      if constexpr ( NDim > 1 ) std::get<1>(_v) = vec[1];
      if constexpr ( NDim > 2 ) std::get<2>(_v) = vec[2];
      return *this;
    }

    constexpr vVec& operator= ( const array<T,N>& arr ) noexcept {
      std::get<0>(_v) = arr[0];
      if constexpr ( NDim > 1 ) std::get<1>(_v) = arr[1];
      if constexpr ( NDim > 2 ) std::get<2>(_v) = arr[2];
      return *this;
    }

  };

}

#endif
