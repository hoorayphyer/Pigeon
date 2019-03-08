#ifndef  _KNL_GRID_HPP_
#define  _KNL_GRID_HPP_

// NOTE convention: the zero of all indices exposed to users start at the first cell in BULK. In other words, guard cells have indices outside [0, dim_of_bulk). Guard is important only when converting from vectorial to linear index, which can be encapsulated in a dedicated function

namespace knl {
  template < typename T >
  struct gridline {
  private:
    T _lower;
    // int _guard;
    int _dim; // doesn't include guards
    T _delta;

  public:
    constexpr gridline() noexcept: gridline( 0.0, 1.0, 0, 1 ) {}

    constexpr gridline( T lower, T upper,
                        // int guard,
                        int dim ) noexcept
      : _lower(lower),
        // _guard(guard),
        _dim(dim),
        _delta( (upper - lower) / dim ) {}

    // abscissa
    constexpr T absc( int i, T shift_from_lb ) const noexcept {
      return _lower + _delta * ( i + shift_from_lb );
    }

    constexpr T lower() const noexcept { return _lower; }
    constexpr T upper() const noexcept { return absc(dim, 0.0); }
    // constexpr int guard() noexcept { return _guard; }
    constexpr int dim() const noexcept { return _dim; }
    constexpr T delta() const noexcept { return _delta; }

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

namespace knl {
  template < typename T >
  struct gridline_slice {
  private:
    const gridline<T>& _gl;
    int _anchor; // the cell in the super gridline, guard cells not included
    int _extent;

  public:
    constexpr gridline_slice( const gridline<T>& gl, int anchor, int extent ) noexcept
      : _gl(gl), _anchor(anchor), _extent(extent) {}
    constexpr T lower() const noexcept { return _gl.absc( _anchor, 0 ); }
    constexpr T upper() const noexcept { return _gl.absc( _anchor + _extent, 0 ); }
    constexpr int dim() const noexcept { return _extent; }
    // constexpr int guard() noexcept { return _gl.guard(); }
    constexpr T delta() const noexcept { return _gl.delta(); }

    constexpr T absc( int i, T shift_from_lb ) const noexcept {
      return _gl( i + _anchor, shift_from_lb );
    }

    constexpr void set_anchor( int anchor ) noexcept { _anchor = anchor; }
    constexpr void set_extent( int extent ) noexcept { _extent = extent; }

  };
}

#include <array>

namespace knl {
  // TODO
  // const std::array< int, Dim_Grid + 1 > stride;

  // template < typename T, int DGrid, template < typename > class GL_t = gridline >
  // using Grid = std::array< GL_t<T>, DGrid >;


  template < typename T, int DGrid, template < typename > class GL_t = gridline >
  struct Grid : public std::array< GL_t<T>, DGrid > {
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
    static constexpr int size = DGrid;

    constexpr Grid( const GL_t<T>& gl0, const GL_t<T>& gl1 ) noexcept
      : std::array< GL_t<T>, DGrid >{ gl0, gl1 } {}

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
  template < int I, typename T, int DGrid, template < typename > class GL_t>
  constexpr const auto& get( const knl::Grid<T,DGrid,GL_t>& grid ) noexcept {
    return std::get<I>( static_cast<const std::array< GL_t<T>, DGrid >&>(grid) );
  }
}

#endif
