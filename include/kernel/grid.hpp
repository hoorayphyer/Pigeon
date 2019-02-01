#ifndef  _KNL_GRID_HPP_
#define  _KNL_GRID_HPP_

#include <array>

namespace knl {
  template < typename T >
  struct gridline_t {
    T lower;
    T upper;
    int guard;
    T delta;

    static_assert( std::is_floating_point_v<T> );

    constexpr gridline_t() noexcept: gridline_t( 0.0, 1.0, 0, 1 ) {}

    constexpr gridline_t( T grid_lower, T grid_upper, int grid_guard, int grid_N ) noexcept
      : lower(grid_lower), upper(grid_upper), guard(grid_guard),
        delta( (grid_upper - grid_lower) / grid_N ) {}

    // abscissa
    constexpr T absc( int i, T shift_from_lb ) noexcept {
      return lower + delta * ( i + shift_from_lb - guard );
    }

    // constexpr T rel_absc( int i, T shift_from_lb ) noexcept {
    //   return (lower / delta) - guard + i + shift_from_lb;
    // }

    // constexpr int index ( T abscissa ) noexcept {
    //   return static_cast<int>( ( abscissa - lower ) / delta + guard );
    // }

    // constexpr int index_from_rel ( T rel_abscissa ) noexcept {
    //   return static_cast<int>( rel_abscissa - (lower / delta) + guard );
    // }
  };

  // Grid is designed to represent the supergrid
  template < std::size_t DGrid, typename T >
  using grid_t = std::array< gridline_t<T>, DGrid >;
}


// TODO replace this with expression template interface
namespace mem {
  template < std::size_t DGrid, typename T >
  auto delta( const knl::grid_t<DGrid,T>& grid ) {
    return std::array<T,DGrid>();
  }

  template < std::size_t DGrid, typename T >
  auto lower( const knl::grid_t<DGrid,T>& grid ) {
    return std::array<T,DGrid>();
  }

  template < std::size_t DGrid, typename T >
  auto upper( const knl::grid_t<DGrid,T>& grid ) {
    return std::array<T,DGrid>();
  }

  template < std::size_t DGrid, typename T >
  auto guard( const knl::grid_t<DGrid,T>& grid ) {
    return std::array<T,DGrid>();
  }
}


#endif
