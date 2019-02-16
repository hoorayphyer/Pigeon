#ifndef  _KNL_GRID_HPP_
#define  _KNL_GRID_HPP_


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

}

#include <array>
#include "apt/virtual_vec.hpp"
#include "apt/vec_from_function.hpp"

namespace knl {
  // Grid is designed to represent the supergrid
  template < int DGrid, typename T >
  struct Grid : public std::array< gridline_t<T>, DGrid >{
    constexpr apt::vVec<T,DGrid> delta() noexcept{
      return apt::make_vff<DGrid>( []( auto& x ) {return x.delta;}, *this);
    }

    constexpr apt::vVec<T,DGrid> lower() noexcept {
      return apt::make_vff<DGrid>( []( auto& x ) {return x.lower;}, *this);
    }

    constexpr apt::vVec<T,DGrid> upper() noexcept {
      return apt::make_vff<DGrid>( []( auto& x ) {return x.upper;}, *this);
    }

    constexpr apt::vVec<int,DGrid> guard() noexcept {
      return apt::make_vff<DGrid>( []( auto& x ) {return x.guard;}, *this);
    }

  };
}


#endif
