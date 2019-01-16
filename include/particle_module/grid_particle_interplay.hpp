#ifndef _GRID_PARTICLE_INTERPLAY_HPP_
#define _GRID_PARTICLE_INTERPLAY_HPP_

#include "core/grid.hpp"
#include "core/vector.hpp"
// #include "field_module/grid_field_interplay.hpp"

namespace particle {
  template < std::size_t DGrid, class Field, class Vec, typename ShapeF >
  auto interpolate_field ( const ShapeF& shape_f, const Field& field, const Grid<DGrid>& grid, Vec location ) {
    static_assert( DGrid < 4 );
    static_assert ( vec::size_v<Vec> >= DGrid );

    // turns location to relative
    location /= mem::delta(grid);

    auto interp_comp =
      [ &location, &grid, &anchor=field.anchor, &shape_f ]( const auto& field_comp ) {
        auto loc_comp = location - field_comp.offset
          - mem::lower(grid) + mem::guard(grid);
        // lower bound index of contributing cells
        auto I_b = vec::per_dim::make<DGrid>
          ( [ r = shape_f.radius ]( const auto& loc ) {
              return int(loc - r) + 1;
            }, loc_comp );

        // separation at I_begin, defined to be loc_comp - I_b, recycle loc_comp
        loc_comp -= I_b;
        const auto& sep_b = loc_comp;

        Real res = 0;
        constexpr auto supp = ShapeF::support();
        // it happens that the one-past-last end = I_b + support
        if constexpr ( DGrid == 1 ) {

            for ( int i = 0; i < supp; ++i ) {
              res += fcomp( i + std::get<0>(I_b) )
                * shape_f( std::get<0>(sep_b) - i );
            }

          } else if ( DGrid == 2 ) {

          for ( int j = 0; j < supp; ++j ) {
            auto wgt_j = shape_f( std::get<1>(sep_b) - j );
            for ( int i = 0; i < supp; ++i ) {
              res += fcomp( i + std::get<0>(I_b),
                            j + std::get<1>(I_b) )
                * shape_f( std::get<0>(sep_b) - i ) * wgt_j;
            }
          }

        } else if ( DGrid == 3 ) {

            for ( int k = 0; k < supp; ++k ) {
              auto wgt_k = shape_f( std::get<2>(sep_b) - k );
              for ( int j = 0; j < supp; ++j ) {
                auto wgt_j = shape_f( std::get<1>(sep_b) - j );
                for ( int i = 0; i < supp; ++i ) {
                  res += fcomp( i + std::get<0>(I_b),
                                j + std::get<1>(I_b),
                                k + std::get<1>(I_b) )
                    * shape_f( std::get<0>(sep_b) - i ) * wgt_j * wgt_k;
                }
              }
            }

        }
        // TODOL
        // else { static_assert( false, "DGrid not implemented" ); }

        return res;
      };

    return vec::per_dim::make<Field::DField>( interp_comp, field.components );
  }

}


#endif
