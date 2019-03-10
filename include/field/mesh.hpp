#ifndef  _FIELD_MESH_HPP_
#define  _FIELD_MESH_HPP_

#include "kernel/grid.hpp"
#include "apt/foreach.hpp"
#include "apt/index.hpp"

namespace field {
  template < typename T, int D >
  struct Mesh {
  private:
    using Bulk_t = knl::Grid<T, D, knl::grid1d::Clip>;

    const Bulk_t& _bulk; // bulk's zero sets the zero of the mesh
    std::array< int[2], D > _margin {}; // margin is where boundary conditions are applied. They don't include guard cells
    int _guard = 0; // uniform value in all directions
    std::array< int, D + 1 > _stride; // bulk and margin combined

  public:
    static constexpr int NDim = D;

    constexpr Mesh( const Bulk_t& bulk, std::array< int[2], NDim > margin, int guard ) noexcept
      : _bulk(bulk), _margin(std::move(margin)), _guard(guard) {
      _stride[0] = 1;
      for ( int i = 0; i < NDim; ++i ) {
        _stride[i+1] = _stride[i] * ( _bulk[i].dim() + _margin[i][0] + _margin[i][1] + 2 * _guard );
      }
    }

    constexpr int linearized_index_of_whole_mesh( const apt::Index<NDim>& i_bulk ) const noexcept {
      // TODO check bounds on i_bulk???
      int I = 0;
      apt::foreach<0,NDim>
        ( [&I,g=_guard]( auto i, auto m, auto s ) {
            I += ( i + m[0] + g ) * s;
          }, i_bulk, _margin, _stride );
      return I;
    }

    constexpr const auto& bulk() const noexcept { return _bulk; }
    constexpr const auto& margin() const noexcept { return _margin; }
    constexpr int guard() const noexcept { return _guard; }

    constexpr int size() const noexcept { return _stride.back(); }

    constexpr apt::Index<NDim> origin() const noexcept { // the first mesh cell expressed in bulk_indices
      apt::Index<NDim> origin;
      apt::foreach<0,NDim>
        ( [g=_guard]( auto& o, const auto& marg ) {
            o = -marg[0] - g;
          }, origin, _margin );
      return origin;
    }

    constexpr apt::Index<NDim> extent() const noexcept { // full extent of the mesh
      apt::Index<NDim> extent;
      for ( int i = 0; i < NDim; ++i ) extent[i] = _stride[i+1] / _stride[i];
      return extent;
    }

  };
}


#endif
