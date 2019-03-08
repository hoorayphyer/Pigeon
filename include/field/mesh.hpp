#ifndef  _FIELD_MESH_HPP_
#define  _FIELD_MESH_HPP_

#include "kernel/grid.hpp"
#include "apt/foreach.hpp"

namespace field {
  template < typename T, int D >
  struct Mesh {
  protected:
    using Bulk_t = knl::Grid<T, D, knl::grid1d::Clip>;

    Bulk_t _bulk; // bulk's zero sets the zero of the manifold
    std::array< int[2], D > _margin;
    std::array< int, D + 1 > _stride; // bulk and margin combined

  public:
    static constexpr int NDim = D;

    constexpr Mesh( const Bulk_t& bulk, const std::array< int[2], NDim >& margin ) noexcept
      : _bulk(bulk), _margin(margin) {
      _stride[0] = 1;
      for ( int i = 0; i < NDim; ++i ) {
        _stride[i+1] = _stride[i] * ( bulk[i].dim() + margin[i][0] + margin[i][1] );
      }
    }

    constexpr int linearized_index_of_whole_mesh( const std::array<int,NDim>& i_bulk ) const noexcept {
      // TODO check bounds on i_bulk???
      int I = 0;
      apt::foreach<0,NDim>
        ( [&I]( auto i, auto m, auto s ) {
            I += ( i + m[0] ) * s;
          }, i_bulk, _margin, _stride );
      return I;
    }

    constexpr int size() const noexcept { return _stride.back(); }
    constexpr const auto& bulk() const noexcept { return _bulk; }
  };
}

#endif
