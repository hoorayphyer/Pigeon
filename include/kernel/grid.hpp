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

namespace knl {
  // Grid here means the supergrid
  template < int DGrid, typename T >
  struct Grid : public std::array< gridline_t<T>, DGrid > {
  private:
    template < typename U >
    using vVec = apt::vVec<const U, DGrid>;

    enum class Mem { DELTA, LOWER, UPPER, GUARD };

    template < Mem M, typename U, std::size_t... I >
    constexpr vVec<U> mem_get( std::index_sequence<I...> ) const noexcept {
      if constexpr ( M == Mem::DELTA )
        return vVec<U>( std::get<I>(*this).delta... );
      else if ( M == Mem::LOWER )
        return vVec<U>( std::get<I>(*this).lower... );
      else if ( M == Mem::UPPER )
        return vVec<U>( std::get<I>(*this).upper... );
      else if ( M == Mem::GUARD )
        return vVec<U>( std::get<I>(*this).guard... );
    }

  public:
    constexpr auto delta() const noexcept {
      return mem_get<Mem::DELTA,T>( std::make_index_sequence<DGrid>{} );
    }

    constexpr auto lower() const noexcept {
      return mem_get<Mem::LOWER,T>( std::make_index_sequence<DGrid>{} );
    }

    constexpr auto upper() const noexcept {
      return mem_get<Mem::UPPER,T>( std::make_index_sequence<DGrid>{} );
    }

    constexpr auto guard() const noexcept {
      return mem_get<Mem::GUARD,int>( std::make_index_sequence<DGrid>{} );
    }

  };
}


#endif
