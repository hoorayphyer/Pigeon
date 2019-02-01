#ifndef _KNL_GRID_SHAPE_INTERPLAY_HPP_
#define _KNL_GRID_SHAPE_INTERPLAY_HPP_

#include "kernel/shape.hpp"
#include "kernel/grid.hpp"
#include "apt/vec.hpp"

namespace knl :: impl {
  template < shape S, typename T, std::size_t DGrid >
  class ShapeRangeInterator {
  private:
    static constexpr shapef_t<S,T> shape_f;

    int _I = 0;

    std::array<T, DGrid - 1> _wgt;
    const Vec<int, DGrid> _I_b;
    const Vec<T, DGrid> _sep_b;

  public:
    using difference_type = int;
    using value_type = void;
    using reference = std::tuple< std::array<int, DGrid>, T >;
    using pointer = void;
    using iterator_category = std::forward_iterator_tag;

    ShapeRangeInterator( int I, const Vec<T,DGrid>& location )
      : _I_b( apt::per_dim::make<DGrid>
              ( []( const auto& loc ) {
                  return int(loc - shape_f.radius) + 1;}, location ) ),
        _sep_b( location - _I_b ) {
    }

    inline bool operator!= ( int I_end ) const {
      return _I != I_end;
    }

    auto& operator++() { ++_I; return *this; }

    // TODO make sure ShapeF can be passed in as constexpr. This may be used to optimize the ever-checking away.
    reference operator*() const {
      constexpr auto supp = sf::support(S);
      int i = _I % supp;
      int j = i / supp;
      int k = i / (supp * supp);
      if ( 0 == i ) {
        wj = shape_f( std::get<1>(sep_b) - j );
        if ( 0 == j ) {
          wk = shape_f( std::get<2>(sep_b) - k );
        }
      }
      return std::make_tuple(_I_b + {i,j,k},
                             shape_f( std::get<0>(sep_b) - i ) * wj * wk );
    }

  };

  template < shape S, typename T, std::size_t DGrid >
  class ShapeRange {
  private:
    Vec<T,DGrid> _loc;
  public:
    ShapeRange( const Vec<T,DGrid>& loc_rel, const grid_t<DGrid>& grid, const Vec<T,DGrid>& offset )
      : _loc ( loc_rel - offset - mem::lower(grid) + mem::guard(grid) ) {}

    auto begin() const && {
                           return ShapeRangeInterator<S, T, DGrid >( 0, std::move(_loc) );
    }

  };

  // TODO double check this function simply returns constexpr. We omitted the argument name
  template < shape S, typename T, std::size_t DGrid >
  constexpr std::end( const ShapeRange< S, T, DGrid > & ) noexcept {
    return DGrid * sf::support(S);
  }
}

namespace knl {
  template < shape S, typename T, std::size_t DGrid >
  inline auto make_shape_range( const Vec<T, DGrid>& loc_rel, const grid_t<DGrid>& grid, const Vec<T,DGrid>& offset ) {
    return ShapeRange<DGrid, S, T>( loc_rel, grid, offset );
  }
}

#endif
