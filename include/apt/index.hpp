#ifndef  _APT_INDEX_HPP_
#define  _APT_INDEX_HPP_

#include "apt/array.hpp"
#include "apt/foreach.hpp"

namespace apt {
  template < int D >
  using Index = array<int,D>;

}

template < int D >
constexpr apt::Index<D> operator+( apt::Index<D> ind_a, const apt::Index<D>& ind_b ) noexcept {
  apt::foreach<0,D>
    ( []( auto& a, auto b ) { a += b; }, ind_a, ind_b );
  return ind_a;
}

#endif
