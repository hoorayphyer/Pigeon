#ifndef  _MAKE_FIELD_HPP_
#define  _MAKE_FIELD_HPP_

#include "field/field.hpp"

namespace field {
  template < typename T, int DField, int DGrid >
  Field< T, DField, DGrid >
  make_field( std::array< int, DGrid > anchor,
              std::array< int, DGrid > extent,
              std::array< std::array< bool, DGrid >, DField > offsets ) {
    Field< T, DField, DGrid > field;
    field.anchor = std::move(anchor);
    field.extent = std::move(extent);

    field.stride[0] = 1;
    for ( int i = 1; i < DGrid + 1; ++i )
      field.stride[i] = field.stride[i-1] * extent[i-1];

    const auto& size = field.stride.back();

    auto itr = offsets.begin();
    for( auto& comp : field ) {
      comp.offset = *(itr++);
      comp.resize(size);
      for ( auto& elm : comp ) elm = static_cast<T>(0);
      comp.shrink_to_fit();
    }
    return field;
  }
}

#endif
