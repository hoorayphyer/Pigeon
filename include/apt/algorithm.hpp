#ifndef  _APT_ALGORITHM_HPP_
#define  _APT_ALGORITHM_HPP_

#include "apt/vec.hpp"

namespace apt :: per_dim {
  namespace impl {

    template < std::size_t D, typename Func, typename... Tuples >
    constexpr auto invoke_on_dim( const Func& f, Tuples&&... tpls ) noexcept {
      return f( std::get<D>(std::forward<Tuples>(tpls))... );
    }

    template < typename Func, typename... Tuples, std::size_t... D >
    constexpr auto make( const Func& f, std::index_sequence<D...>, Tuples&&... tpls ) noexcept {
      return apt::make( invoke_on_dim<D>( f, std::forward<Tuples>(tpls)... )... );
    }

    template < typename Func, typename... Tuples, std::size_t... D >
    constexpr auto tie( const Func& f, std::index_sequence<D...>, Tuples&... tpls ) noexcept {
      return apt::tie( invoke_on_dim<D>( f, tpls... )... );
    }

  }


  // f shouldn't depend on the order of each dim executed
  // TODO return Vec<T> or Vec<T&> ??? How to decide Vec<T>& vs Vec<T&>. Maybe just return Vec<T> or Vec<T&> ???
  // NOTE TODO currently per_dim_of requires Func to return something, and the return object of per_dim_of is always apt::made rather than apt::tied, so there is always a temporary object returned. Cannot use per_dim_of to generate reference
  template < std::size_t Ndims, typename Func, typename... Vectors >
  constexpr auto make( const Func& f, Vectors&&... vecs  ) noexcept {
    return impl::make( f, std::make_index_sequence<Ndims>{},
                       std::forward<Vectors>(vecs)... );
  }

  template < std::size_t Ndims, typename Func, typename... Vectors >
  constexpr auto tie( const Func& f, Vectors&... vecs  ) noexcept {
    return impl::tie( f, std::make_index_sequence<Ndims>{}, vecs... );
  }

}

namespace apt {
  // NOTE foreach has no return type as opposed to those in apt::per_dim
  template < std::size_t Begin, std::size_t End, typename Func, typename... Vectors >
  constexpr void foreach( const Func& f, Vectors&&... vecs  ) noexcept {
    static_assert( Begin <= End );
    if constexpr ( Begin == End ) return;
    else {
      f( std::get<Begin>(std::forward<Vectors>(vecs))... );
      return foreach<Begin+1, End>( f, std::forward<Vectors>(vecs)... );
    }
  }
}

#include <algorithm> // for std::swap
namespace apt {
  template < typename T, std::size_t N >
  void swap( Vec<T,N>& a, Vec<T,N>& b ) noexcept {
    foreach<0,N>
      ( []( auto& x, auto& y ) noexcept { std::swap(x,y); }, a, b );
  }
}

#endif
