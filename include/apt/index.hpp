#ifndef  _APT_INDEX_HPP_
#define  _APT_INDEX_HPP_

namespace apt {
  template < int D >
  struct Index {
    int data[D] {};

    constexpr int operator[] ( int i ) const noexcept { return data[i]; }
    constexpr int& operator[] ( int i ) noexcept { return data[i]; }

    constexpr Index operator+( Index other ) const noexcept;
  };

}

namespace std {
  // define this so as to be used in apt::foreach
  template < int I, int D >
  constexpr auto get( const apt::Index<D>& index ) noexcept {
    static_assert( I < D );
    return index[I];
  }

  template < int I, int D >
  constexpr auto& get( apt::Index<D>& index ) noexcept {
    static_assert( I < D );
    return index[I];
  }
}

#include "apt/foreach.hpp"
namespace apt {
  template < int D >
  constexpr Index<D> Index<D>::operator+( Index<D> other ) const noexcept {
    foreach<0,D>( [](int& a, int b ) { a += b; }, other, *this );
    return other;
  }

}

#endif
