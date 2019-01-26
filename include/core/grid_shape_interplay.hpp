#ifndef _GRID_SHAPE_INTERPLAY_HPP_
#define _GRID_SHAPE_INTERPLAY_HPP_

#include "core/shapefunction.hpp"
#include "core/grid.hpp"
#include "core/vector.hpp"

namespace sf :: impl {
  template < std::size_t DGrid, sf::shape S, typename T >
  class ShapeRangeInterator {
  private:
    static constexpr sf::ShapeFunction<S,T> shape_f;

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

  template < std::size_t DGrid, sf::shape S, typename T >
  class ShapeRange {
  private:
    Vec<T,DGrid> _loc;
  public:
    ShapeRange( const Vec<T,DGrid>& loc_rel, const Grid<DGrid>& grid, const Vec<T,DGrid>& offset )
      : _loc ( loc_rel - offset - mem::lower(grid) + mem::guard(grid) ) {}

    auto begin() const && {
                           return ShapeRangeInterator<DGrid, S, T>( 0, std::move(_loc) );
    }

  };

  // TODO double check this function simply returns constexpr. We omitted the argument name
  template < std::size_t DGrid, sf::shape S, typename T >
  constexpr std::end( const ShapeRange<DGrid, S, T> & ) noexcept {
    return DGrid * sf::support(S);
  }
}

namespace sf {
  template < std::size_t DGrid, sf::shape S, typename T >
  inline auto make_shape_range( const Vec<T, DGrid>& loc_rel, const Grid<DGrid>& grid, const Vec<T,DGrid>& offset ) {
    return ShapeRange<DGrid, S, T>( loc_rel, grid, offset );
  }
}

#endif
