#ifndef  _INTERPOLATION_HPP_
#define  _INTERPOLATION_HPP_

#include "kernel/grid_shape_interplay.hpp"
#include "field/field.hpp"

namespace field {

  template < typename T, int DField, std::size_t DGrid, class Vec, typename ShapeF >
  auto interpolate ( const Field<T,DField,DGrid>& field, const Vec& location,
                     const knl::Grid<DGrid, T>& grid, const ShapeF& shapef ) {
    // only need first DGrid elements of location even if location is of higher dim than DGrid
    // TODO expression template
    auto loc_rel = apt::per_dim::make<DGrid>
      ( []( const auto& l, const auto& delta ) {
          return l / delta;
        }, location, mem::delta(grid) );

    auto interp_comp =
      [ &loc_rel, &grid, &anchor=field.anchor ] ( const auto& f_comp ) {
        T res = 0;
        for ( auto[ I, wgt ] : knl::make_shape_range( loc_rel, grid, f_comp.offset ) )
          res += field_comp( I ) * wgt;
        return res;
      };

    return apt::per_dim::make<DField>( interp_comp, field );
  }

}

#endif
