#ifndef  _KNL_GRID_HPP_
#define  _KNL_GRID_HPP_

#include "kernel/grid_1d.hpp"
#include "apt/array.hpp"
#include "apt/index.hpp"

namespace knl {
  template < typename T, int DGrid, template < typename > class grid1d = grid1d::Whole >
  struct Grid : public apt::array< grid1d<T>, DGrid > {
    using element_type = T;
    static constexpr int NDim = DGrid;

    using apt::array< grid1d<T>, DGrid >::array;

    constexpr apt::Index<DGrid> extent() const noexcept {
      apt::Index<DGrid> ext;
      for ( int i = 0; i < DGrid; ++i )
        ext[i] = (*this)[i].dim();
      return ext;
    }

  };

  template < typename T, int DGrid >
  Grid<T,DGrid,grid1d::Clip> make_clip( const Grid<T,DGrid,grid1d::Whole>& grid ) {
    static_assert( DGrid > 0 && DGrid < 4 );
    // if constexpr ( DGrid == 1 )
    //   return Grid<T,DGrid,grid1d::Clip>(grid[0]);
    // else if ( DGrid == 2 )
    // TODO TODO
    // return Grid<T,DGrid,grid1d::Clip>{{ grid1d::Clip(grid[0]), grid1d::Clip(grid[1]) }};
    // else
    //   return Grid<T,DGrid,grid1d::Clip>(grid[0], grid[1], grid[2]);
  }


}

#endif
