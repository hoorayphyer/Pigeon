#ifndef  _FIELD_MESH_HPP_
#define  _FIELD_MESH_HPP_

#include "kernel/grid.hpp"
#include "apt/index.hpp"
#include "apt/pair.hpp"

namespace field {
  template < typename T, int D >
  struct Mesh {
  private:
    using Bulk_t = knl::Grid<T, D, knl::grid1d::Clip>;

    const Bulk_t& _bulk; // bulk's zero sets the zero of the mesh
    apt::array< apt::pair<int>, D > _margin {}; // margin is where boundary conditions are applied. They don't include guard cells
    int _guard = 0; // uniform value in all directions

  public:
    static constexpr int NDim = D;

    constexpr Mesh( const Bulk_t& grid_bulk, int guard ) noexcept
      : _bulk(grid_bulk), _guard(guard) {}

    constexpr int linearized_index_of_whole_mesh( const apt::Index<NDim>& i_bulk ) const noexcept {
      // TODO check bounds on i_bulk???
      int i = NDim - 1;
      int I = i_bulk[i] + _margin[i][0] + _guard;
      for ( --i; i > -1; --i )
        I = ( i_bulk[i] + _margin[i][0] + _guard ) + I * extent(i);
      return I;
    }

    constexpr const auto& bulk() const noexcept { return _bulk; }
    constexpr const auto& margin() const noexcept { return _margin; }
    constexpr int guard() const noexcept { return _guard; }

    constexpr apt::Index<NDim> origin() const noexcept { // the first mesh cell expressed in bulk_indices
      apt::Index<NDim> origin;
      for ( int i = 0; i < NDim; i++ )
        origin[i] = - _margin[i][0] - _guard;
      return origin;
    }

    constexpr int extent( int ith_dim ) const noexcept { // full extent of the mesh
      return _bulk[ith_dim].dim() + _margin[ith_dim][0] + _margin[ith_dim][1] + 2 * _guard;
    }

    constexpr apt::Index<NDim> extent() const noexcept { // full extent of the mesh
      apt::Index<NDim> ext;
      for ( int i = 0; i < NDim; ++i )
        ext[i] = extent(i);
      return ext;
    }

  };
}


#endif
