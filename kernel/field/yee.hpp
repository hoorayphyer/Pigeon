#ifndef  _FIELD_YEE_HPP_
#define  _FIELD_YEE_HPP_

#include "field/offset.hpp"

namespace field::yee {
  constexpr offset_t Etype = INSITU; // i.e. E[i] is INSITU in i-th dim
  constexpr offset_t Btype = MIDWAY;

  template < int DGrid >
  constexpr auto ofs_gen( offset_t type, int comp ) noexcept {
    apt::array<field::offset_t, DGrid> res;
    apt::foreach<0,DGrid>( [&type]( auto& ofs ) { ofs = {!type}; }, res );
    if ( comp < DGrid ) res[comp] = {type};
    return res;
  }

  template < offset_t Ftype >
  constexpr offset_t ofs_get( int comp, int grid_dim ) noexcept { // NOTE: order doesn't matter
    return ( comp == grid_dim ) ? Ftype : !Ftype;
  }
}

#endif
