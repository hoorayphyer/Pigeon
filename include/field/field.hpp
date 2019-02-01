#ifndef  _FIELD_HPP_
#define  _FIELD_HPP_

#include "apt/numeric.hpp"
#include <vector>

namespace field {
  template < typename T, int DGrid >
  struct Component : std::vector<T> {
    const std::array< T, DGrid > offset;
  };
}

namespace field {
  template < typename T, int Dim_Field, int Dim_Grid >
  struct Field : public std::array < Component<T, Dim_Grid>, Dim_Field > {
    using data_type = T;
    static constexpr auto DField = Dim_Field;
    static constexpr auto DGrid = Dim_Grid;

    const std::array< int, Dim_Grid > anchor;
    const std::array< int, Dim_Grid > stride;

    // //------ int... inds_local -------
    // // TODO indices
    // template < int Comp >
    // inline decltype(auto) c( int... inds ) {
    //   static_assert( sizeof...(inds) >= Field::DGrid, "Not enough indices provided" );
    //   // TODO inner_product
    //   auto&& ind_linear = tum::inner_product( field.stride, std::forward_as_tuple(std::move(inds)...) );
    //   return std::get<Comp>(components)[std::move(ind_linear)];
    // }

    // template < int Comp >
    // inline decltype(auto) c( int... inds ) const {
    //   static_assert( sizeof...(inds) >= Field::DGrid, "Not enough indices provided" );
    //   // TODO inner_product
    //   auto&& ind_linear = tum::inner_product( field.stride, std::forward_as_tuple(std::move(inds)...) );
    //   return std::get<Comp>(components)[std::move(ind_linear)];
    // }

    //------ int... inds_global -------
    template < int Comp >
    inline decltype(auto) c( const std::array<int,DGrid>& I_global ) {
      // TODO check if apt::dot can work on std::array
      auto&& ind_linear = apt::dot( stride, I_global - anchor );
      return std::get<Comp>(*this)[std::move(ind_linear)];
    }

    template < int Comp >
    inline const decltype(auto) c( const std::array<int,DGrid>& I_global ) const {
      auto&& ind_linear = apt::dot( stride, I_global - anchor );
      return std::get<Comp>(*this)[std::move(ind_linear)];
    }

  };
}

namespace field {
  template < typename T, int DField, int DGrid >
  Field< T, DField, DGrid >
  make_field( std::array< int, DGrid > anchor,
              std::array< int, DGrid > extent,
              std::array< std::array< bool, DGrid >, DField > offsets );
}

#endif
