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
  auto get( Field&& field, Ind... indices ) {
    static_assert( sizeof...(indices) >= Field::DGrid, "Not enough indices provided" );
    // TODOL check if the second array which is tmp can be used here
    return tum::map( [li=tum::inner_product( field.stride, std::make_tuple(indices...) )]
                     ( auto&& comp ){return comp[li];}, field.comp );
  }





}


// TODO fix this when doing field updater
// void operator*= ( const T value ) {
//   auto f = []( auto& elm ) { elm *= value; };
//   if constexpr ( Field_Dim == 1 ) {
//     std::for_each( _data[0].begin(), _data[0].end(), f );
//   }
//   else if ( Field_Dim == 3 ) {
//     std::for_each( _data[0].begin(), _data[0].end(), f );
//     std::for_each( _data[1].begin(), _data[1].end(), f );
//     std::for_each( _data[2].begin(), _data[2].end(), f );
//   }
// }

// void operator+= ( const Field& field ) {}

#Endif
