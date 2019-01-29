#ifndef _GRID_PARTICLE_INTERPLAY_HPP_
#define _GRID_PARTICLE_INTERPLAY_HPP_

#include "core/grid_shape_interplay.hpp"

namespace particle {

  template < std::size_t DGrid,  class Field, class Vec >
  auto interpolate_field ( const Field& field, const Grid<DGrid>& grid, const Vec& location ) {
    // only need first DGrid elements of location even if location is of higher dim than DGrid
    auto loc_rel = apt::per_dim::make<DGrid>
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

    return apt::per_dim::make<Field::DField>( interp_comp, field.components );
  }

}

#endif
