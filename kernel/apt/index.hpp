#ifndef  _APT_INDEX_HPP_
#define  _APT_INDEX_HPP_

#include "apt/array.hpp"
#include "apt/foreach.hpp"

namespace apt {
  template < int D >
  using Index = array<int,D>;

  template < int D >
  constexpr auto trans(int dim, Index<D> Ib, Index<D> extent ) noexcept {
    Ib[dim] = 0;
    extent[dim] = 1;
    return std::make_pair(Ib,extent);
  }
}

template < int D >
constexpr apt::Index<D> operator+( const apt::Index<D>& ind_a, const apt::Index<D>& ind_b ) noexcept {
  apt::Index<D> res;
  apt::foreach<0,D>
    ( []( auto& r, auto a, auto b ) { r = a + b; }, res, ind_a, ind_b );
  return res;
}

template < int D >
constexpr apt::Index<D> operator-( const apt::Index<D>& ind_a, const apt::Index<D>& ind_b ) noexcept {
  apt::Index<D> res;
  apt::foreach<0,D>
    ( []( auto& r, auto a, auto b ) { r = a - b; }, res, ind_a, ind_b );
  return res;
}

template < int D >
constexpr apt::Index<D>& operator+=( apt::Index<D>& ind_a, const apt::Index<D>& ind_b ) noexcept {
  apt::foreach<0,D>
    ( []( auto& a, auto b ) { a += b; }, ind_a, ind_b );
  return ind_a;
}

template < int D >
constexpr bool operator==( const apt::Index<D>& ind_a, const apt::Index<D>& ind_b ) noexcept {
  bool res = true;
  apt::foreach<0,D>
    ( [&res]( auto a, auto b ) { res = ( res && (a == b) ); }, ind_a, ind_b );
  return res;
}

template < int D >
constexpr bool operator!=( const apt::Index<D>& ind_a, const apt::Index<D>& ind_b ) noexcept {
  return !( ind_a == ind_b );
}

#endif
