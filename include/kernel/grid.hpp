#ifndef  _KNL_GRID_HPP_
#define  _KNL_GRID_HPP_

#include "kernel/grid_1d.hpp"
#include "apt/array.hpp"

namespace knl {
  template < typename T, int DGrid, template < typename > class grid1d = grid1d::Whole >
  struct Grid : public apt::array< grid1d<T>, DGrid > {
    using element_type = T;
    static constexpr int NDim = DGrid;

    // TODOH use DGrid as number of arguments. Or use factory?
    constexpr Grid( const grid1d<T>& gl0, const grid1d<T>& gl1 ) noexcept
      : apt::array< grid1d<T>, DGrid >{ gl0, gl1 } {}

    constexpr apt::Index<DGrid> extent() const noexcept {
      apt::Index<DGrid> ext;
      for ( int i = 0; i < DGrid; ++i )
        ext[i] = (*this)[i].dim();
      return ext;
    }

  };
}

#endif
