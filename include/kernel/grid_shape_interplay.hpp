#ifndef _KNL_GRID_SHAPE_INTERPLAY_HPP_
#define _KNL_GRID_SHAPE_INTERPLAY_HPP_

#include "kernel/shape.hpp"
#include "kernel/grid.hpp"
#include "apt/vec.hpp"

namespace knl :: impl {
  template < typename T, std::size_t DGrid, typename ShapeF >
  class ShapeRangeInterator {
  private:
    const ShapeF& _sf;

    int _I = 0;

    std::array<int, DGrid> _ijk;
    std::array<T, DGrid> _wgt;
    const apt::Vec<int, DGrid> _I_b;
    const apt::Vec<T, DGrid> _sep_b;

  public:
    using difference_type = int;
    using value_type = void;
    using reference = std::tuple< std::array<int, DGrid>, T >;
    using pointer = void;
    using iterator_category = std::forward_iterator_tag;

    ShapeRangeInterator( int I, const apt::Vec<T,DGrid>& location, const ShapeF& shapef )
      : _sf(shapef),
        _I_b( apt::make_vff<DGrid>
              ( [&sf=shapef]( const auto& loc ) {
                  return int(loc - sf.support / 2.0) + 1;}, location ) ),
        _sep_b( _I_b - location ) {}

    inline bool operator!= ( int I_end ) const {
      return _I != I_end;
    }

    auto& operator++() { ++_I; return *this; }

    // TODO make sure ShapeF can be passed in as constexpr. This may be used to optimize the ever-checking away.
    reference operator*() const {
      std::get<0>(_ijk) = _I % _sf.support;
      std::get<0>(_wgt) = _sf( std::get<0>(_ijk) + std::get<0>(_sep_b) );

      if constexpr ( DGrid >= 2 ) {
          if( std::get<0>(_ijk) != 0 ) goto W1; // NOTE label should be unique inside a function scope
          std::get<1>(_ijk) = ( _I % ( _sf.support * _sf.support ) ) / _sf.support;
          std::get<1>(_wgt) = _sf( std::get<1>(_ijk) + std::get<1>(_sep_b) );
        W1:
          std::get<0>(_wgt) *= std::get<1>(_wgt);
          goto FINISH;
      }

      if constexpr ( DGrid == 3 ) {
          if( std::get<1>(_ijk) != 0 ) goto W2;
          std::get<2>(_ijk) = _I / ( _sf.support * _sf.support );
          std::get<2>(_wgt) = _sf( std::get<2>(_ijk) + std::get<2>(_sep_b) );
        W2:
          std::get<0>(_wgt) *= std::get<2>(_wgt);
          goto FINISH;
      }

    FINISH:
      return std::make_tuple( _I_b + _ijk, std::get<0>(_wgt) );
    }

  };

  template < typename T, std::size_t DGrid, typename ShapeF >
  class ShapeRange {
  private:
    apt::Vec<T,DGrid> _loc;
    const ShapeF& _sf;

  public:
    ShapeRange( const apt::Vec<T,DGrid>& loc_rel, const Grid<DGrid, T>& grid,
                const apt::Vec<T,DGrid>& offset, const ShapeF& shapef )
      : _loc ( loc_rel - offset - grid.lower() + grid.guard() ), _sf(shapef) {}

    auto begin() const && {
                           return ShapeRangeInterator<T, DGrid, ShapeF>( 0, std::move(_loc), _sf );
    }

    auto end() const && noexcept { return DGrid * _sf.support; }
  };

}

namespace knl {
  template < typename T, std::size_t DGrid, typename ShapeF >
  inline auto make_shape_range( const apt::Vec<T, DGrid>& loc_rel, const Grid<DGrid, T>& grid,
                                const apt::Vec<T,DGrid>& offset, const ShapeF& shapef ) {
    return impl::ShapeRange<T, DGrid, ShapeF>( loc_rel, grid, offset, shapef );
  }
}

#endif
