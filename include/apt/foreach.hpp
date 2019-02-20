#ifndef _APT_FOREACH_HPP_
#define _APT_FOREACH_HPP_

#include <utility> // for std::forward

// declare std::get for std::array
namespace std {
  template < typename, std::size_t > struct array;
  template< size_t I, class T, size_t N >
  constexpr T& get( array<T,N>& a ) noexcept;

  template< size_t I, class T, size_t N >
  constexpr T&& get( array<T,N>&& a ) noexcept;

  template< size_t I, class T, size_t N >
  constexpr const T& get( const array<T,N>& a ) noexcept;

  template< size_t I, class T, size_t N >
  constexpr const T&& get( const array<T,N>&& a ) noexcept;
}

namespace apt {
  template < std::size_t Begin, std::size_t End, typename Func, typename... Args >
  constexpr void foreach( const Func& f, Args&&... args  ) noexcept {
    static_assert( Begin <= End );
    if constexpr ( Begin == End ) return;
    else {
      f( std::get<Begin>(std::forward<Args>(args))... );
      return foreach<Begin+1, End>( f, std::forward<Args>(args)... );
    }
  }

}

#endif
