#ifndef  _FIELD_HPP_
#define  _FIELD_HPP_

#include "phys_vector.hpp"
#include "grid.hpp"
#include <array>
#include <vector>
#include <tuple>

// Field indexing convention: ( 0, 0, 0 ) is the origin of the bulk.

// #include <algorithm> // std::for_each, std::transform

// Manifold_Dim is here more for distinction from Field_Dim.
// TODOL look at fold-expression in C++17, and tuple index_sequences
template < typename T, int Field_Dim, int Manifold_Dim, bool Midpoint1, bool Midpoint2, bool Midpoint3 >
class Field {
private:
  std::array< std::vector<T>, Field_Dim > _data;

  constexpr int linear_index( const int i, const int j, const int k ) const {
    return i + j * grid_patch<0>::extent_with_guard() + k * grid_patch<0>::extent_with_guard() * grid_patch<1>::extent_with_guard() +
      grid<0>::guard + grid<1>::guard * grid_patch<0>::extent_with_guard() + grid<2>::guard * grid_patch<0>::extent_with_guard() * grid_patch<1>::extent_with_guard();
  }

public:
  Field( const std::array<int, Manifold_Dim> extent_with_guard ) {
    static_assert( Manifold_Dim == 3, "Sorry, only 3D is implemented for now");
    if constexpr( Manifold_Dim == 3 ) {
      int size = extent_with_guard[0] * extent_with_guard[1] * extent_with_guard[2];
      for ( auto& elm : _data ) {
        elm = std::vector<T>( size, T() );
        elm.shrink_to_fit();
      }
    }
  }

  Field() : Field( { grid_patch<0>::extent_with_guard(), grid_patch<1>::extent_with_guard(), grid_patch<2>::extent_with_guard() } ) {}

  Field( const Field& other ) = default;
  Field( Field&& other ) = default;

  // NOTE decltype(auto) is necessary
  inline decltype(auto) operator() ( int i, int j = 0, int k = 0 ) {
    if constexpr ( Field_Dim == 1 ) return _data[ linear_index(i, j, k) ];
    else if ( Field_Dim == 3 ) {
      i = linear_index( i, j, k );
      return std::forward_as_tuple( _data[0](i), _data[1](i), _data[2](i)  );
    }
  }

  inline decltype(auto) operator() ( int i, int j = 0, int k = 0 ) const {
    if constexpr ( Field_Dim == 1 ) return _data[ linear_index(i, j, k) ];
    else if ( Field_Dim == 3 ) {
      i = linear_index( i, j, k );
      return std::forward_as_tuple( _data[0](i), _data[1](i), _data[2](i)  );
    }
  }

  inline auto atlas( int i, int j = 0, int k = 0 ) const {
    return std::make_tuple(
                           grid<0>::abscissa<Midpoint1>( grid_patch<0>::origin + i ),
                           grid<1>::abscissa<Midpoint2>( grid_patch<1>::origin + j ),
                           grid<2>::abscissa<Midpoint3>( grid_patch<2>::origin + k )
                           );

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

};




#endif
