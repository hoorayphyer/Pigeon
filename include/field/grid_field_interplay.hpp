#ifndef _GRID_FIELD_INTERPLAY_HPP_
#define _GRID_FIELD_INTERPLAY_HPP_

#include "core/field.hpp"
#include "core/grid.hpp"
#include "utility/tuple_manip.hpp"
#include <numeric>
#include <algorithm>

namespace field {
  template < typename T, int DField, int DGrid >
  auto make_field( std::array< int, DGrid > anchor, std::array< int, DGrid > extent, std::array< std::array< bool, DGrid > > offsets ) {
    auto stride = extent;
    std::exclusive_scan(extent.begin(), extent.end(), stride.begin(), 1, [](int a, int b){return a*b; } );
    auto size = stride.back() * extent.back();
    return Field<T, DField, DGrid>
      { anchor, stride,
          tum::map( [size=size] (auto offset)
                    { return FieldComponent<T,DGrid>(size, offset)}, offsets ) };
  }

  template <int Comp=0, class Field, int... Ind>
  inline auto val( Field&& field, Ind... indices ) {
    static_assert( sizeof...(indices) >= Field::DGrid, "Not enough indices provided" );
    auto&& ind_linear = tum::inner_product( field.stride, std::forward_as_tuple(indices...) );
    return field[Comp][std::move(ind_linear)];
  }

  template <int Comp=0, class Field, int... Ind>
  inline auto atlas( Field&& field, const std::array<Grid, Field::DGrid>& manifold, Ind... indices ) {
    static_assert( sizeof...(indices) >= Field::DGrid, "Not enough indices provided" );
    return tum::map([](auto&& grid, auto&& index, auto&& offset, auto&& anchor) {
                      return grid.abscissa( index, offset + anchor ) },
                    manifold, std::forward_as_tuple(indices), field[Comp].offset, field.anchor );
  }



}

#Endif
