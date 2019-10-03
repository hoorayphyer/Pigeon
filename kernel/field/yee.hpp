#ifndef  _FIELD_YEE_HPP_
#define  _FIELD_YEE_HPP_

#include "apt/array.hpp"
#include "field/offset.hpp"

namespace field::yee {
  constexpr offset_t Etype = INSITU; // i.e. E[i] is INSITU in i-th dim
  constexpr offset_t Btype = MIDWAY;

  template < int DGrid >
  constexpr auto ofs_gen( offset_t type, int comp ) noexcept {
    apt::array<field::offset_t, DGrid> res;
    for ( int i = 0; i < DGrid; ++i ) res[i] = {!type};
    if ( comp < DGrid ) res[comp] = {type};
    return res;
  }
}

#endif
