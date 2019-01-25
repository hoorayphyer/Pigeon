#include "core/field.hpp"

template < typename T, int DField, int DGrid >
Field< T, DField, DGrid >
make_field( std::array< int, DGrid > anchor, std::array< int, DGrid > extent, std::array< std::array< bool, DGrid >, DField > offsets ) {

  auto stride = extent;
  std::exclusive_scan(extent.begin(), extent.end(), stride.begin(), 1, [](int a, int b){return a*b; } );
  auto size = stride.back() * extent.back();

  // TODOL per_dim on std::array missing
  if constexpr ( DField == 1 ) {
    return { std::move(anchor), std::move(stride),
             { FieldComponent<T,DGrid>(size, std::get<0>(offsets) ) } };
  } else if ( DField == 2 ) {
    return { std::move(anchor), std::move(stride),
             { FieldComponent<T,DGrid>(size, std::get<0>(offsets) ),
               FieldComponent<T,DGrid>(size, std::get<1>(offsets) ) } };
  } else if ( DField == 3 ) {
    return { std::move(anchor), std::move(stride),
             { FieldComponent<T,DGrid>(size, std::get<0>(offsets) ),
               FieldComponent<T,DGrid>(size, std::get<1>(offsets) ),
               FieldComponent<T,DGrid>(size, std::get<2>(offsets) ) } };
  }
}
