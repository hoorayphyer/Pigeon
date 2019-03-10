#ifndef  _KNL_GRID_HPP_
#define  _KNL_GRID_HPP_

#include "kernel/grid_1d.hpp"
#include <array>

namespace knl {
  template < typename T, int DGrid, template < typename > class grid1d = grid1d::Whole >
  struct Grid : public std::array< grid1d<T>, DGrid > {
  private:
    static_assert( std::is_floating_point_v<T> );

    // enum class Mem { DELTA, LOWER, UPPER, DIM };

    // template < Mem M, typename U, std::size_t... I >
    // constexpr std::array<U,DGrid> mem_get( std::index_sequence<I...> ) const noexcept {
    //   if constexpr ( M == Mem::DELTA )
    //     return { std::get<I>(*this).delta()... };
    //   else if ( M == Mem::LOWER )
    //     return { std::get<I>(*this).lower()... };
    //   else if ( M == Mem::UPPER )
    //     return { std::get<I>(*this).upper()... };
    //   else if ( M == Mem::DIM )
    //     return { std::get<I>(*this).dim()... };
    // }

  public:
    using element_type = T;
    static constexpr int NDim = DGrid;

    // TODO use DGrid as number of arguments
    constexpr Grid( const grid1d<T>& gl0, const grid1d<T>& gl1 ) noexcept
      : std::array< grid1d<T>, DGrid >{ gl0, gl1 } {}

    constexpr const auto& operator[] ( int i ) const noexcept {
      return std::array< grid1d<T>, DGrid >::operator[] (i);
    }

    // constexpr Grid( const GL_t<T>& gl0, const GL_t<T>& gl1 ) noexcept
    //   : std::array< GL_t<T>, DGrid >{ gl0, gl1 } {}

    // constexpr auto deltas() const noexcept {
    //   return mem_get<Mem::DELTA,T>( std::make_index_sequence<DGrid>{} );
    // }

    // constexpr auto lowers() const noexcept {
    //   return mem_get<Mem::LOWER,T>( std::make_index_sequence<DGrid>{} );
    // }

    // constexpr auto uppers() const noexcept {
    //   return mem_get<Mem::UPPER,T>( std::make_index_sequence<DGrid>{} );
    // }

    // constexpr auto dims() const noexcept {
    //   return mem_get<Mem::DIM,int>( std::make_index_sequence<DGrid>{} );
    // }

  };
}

namespace std {
  // define this so as to be used in apt::foreach
  template < int I, typename T, int DGrid, template < typename > class G>
  constexpr const auto& get( const knl::Grid<T,DGrid,G>& grid ) noexcept {
    return std::get<I>( static_cast<const std::array< G<T>, DGrid >&>(grid) );
  }
}

#endif
