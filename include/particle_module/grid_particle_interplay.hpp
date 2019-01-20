#ifndef _GRID_PARTICLE_INTERPLAY_HPP_
#define _GRID_PARTICLE_INTERPLAY_HPP_

#include "core/grid_shape_interplay.hpp"
// #include "field_module/grid_field_interplay.hpp"

namespace particle {

  template < std::size_t DGrid,  class Field, class Vec >
  auto interpolate_field ( const Field& field, const Grid<DGrid>& grid, const Vec& location ) {
    // only need first DGrid elements of location even if location is of higher dim than DGrid
    auto loc_rel = vec::per_dim::make<DGrid>
      ( []( const auto& l, const auto& delta ) {
          return l / delta;
        }, location, mem::delta(grid) );

    auto interp_comp =
      [ &loc_rel, &grid, &anchor=field.anchor, &shape_f ] ( const auto& f_comp ) {
        Real res = 0;
        for ( auto[ I, wgt ] : sf::make_shape_range( loc_rel, grid, f_comp.offset ) )
          res += field_comp( I ) * wgt;
        return res;
      };

    return vec::per_dim::make<Field::DField>( interp_comp, field.components );
  }

}


// // TODO use iterator as generator to return a range
// auto shape_range = // TODO check auto location, I want a copy all the time
//   []( auto location, const auto& grid, const auto& offset, const auto& shape_f ) {
//     // turns location to relative
//     location /= mem::delta(grid);
//     // TODO check why is there no anchor, is it because
//     location -= ( offset + mem::lower(grid) - mem::guard(grid) );
//     // lower bound index of contributing cells
//     auto I_b = vec::per_dim::make<DGrid>
//       ( [ r = shape_f.radius ]( const auto& loc ) {
//           return int(loc - r) + 1;
//         }, location );

//     location -= I_b; // recycle location to mean separation at I_begin, defined to be location - I_b
//     return std::make_tuple( std::move(I_b), std::move(location) );
//   };

  // constexpr auto supp = ShapeF::support();
  // it happens that the one-past-last end = I_b + support
  //   if constexpr ( DGrid == 1 ) {

  //       for ( int i = 0; i < supp; ++i ) {
  //         res += fcomp( i + std::get<0>(I_b) )
  //           * shape_f( std::get<0>(sep_b) - i );
  //       }

  //     } else if ( DGrid == 2 ) {

  //     for ( int j = 0; j < supp; ++j ) {
  //       auto wgt_j = shape_f( std::get<1>(sep_b) - j );
  //       for ( int i = 0; i < supp; ++i ) {
  //         res += fcomp( i + std::get<0>(I_b),
  //                       j + std::get<1>(I_b) )
  //           * shape_f( std::get<0>(sep_b) - i ) * wgt_j;
  //       }
  //     }

  //   } else if ( DGrid == 3 ) {

  //       for ( int k = 0; k < supp; ++k ) {
  //         auto wgt_k = shape_f( std::get<2>(sep_b) - k );
  //         for ( int j = 0; j < supp; ++j ) {
  //           auto wgt_j = shape_f( std::get<1>(sep_b) - j );
  //           for ( int i = 0; i < supp; ++i ) {
  //             res += fcomp( i + std::get<0>(I_b),
  //                           j + std::get<1>(I_b),
  //                           k + std::get<1>(I_b) )
  //               * shape_f( std::get<0>(sep_b) - i ) * wgt_j * wgt_k;
  //           }
  //         }
  //       }

  //   } else { static_assert( DGrid < 4, "not implemented" ); }


#endif
