#ifndef  _KNL_GRID_HPP_
#define  _KNL_GRID_HPP_

// NOTE convention: the zero of all indices exposed to users start at the first cell in BULK. In other words, guard cells have indices outside [0, dim_of_bulk). Guard is important only when converting from vectorial to linear index, which can be encapsulated in a dedicated function

// guard, indent are controlled by fields directly, where they are collectively called margin cells.

namespace knl :: grid1d {
  template < class E >
  struct Expression {
  private:
    constexpr auto _lower_rel() const noexcept { return static_cast<const E&>(*this)._lower_rel(); }

  public:
    constexpr int dim() const noexcept { return static_cast<const E&>(*this).dim(); }
    constexpr auto delta() const noexcept { return static_cast<const E&>(*this).delta(); }

    constexpr auto lower() const noexcept { return _lower_rel() * delta(); }
    constexpr auto upper() const noexcept { return absc(dim()); }

    // abscissa
    template < typename U >
    constexpr auto absc( int i, U shift_from_lb = 0.0 ) const noexcept {
      return  delta() * ( _lower_rel() + i + shift_from_lb );
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

namespace knl :: grid1d {
  template < typename T > struct Clip;

  template < typename T >
  struct Whole : Expression<Whole<T>> {
  private:
    int _dim;
    T _delta;
    T _lo_rel;

    constexpr T _lower_rel() const noexcept { return _lo_rel; }
  public:
    friend class Expression<Whole<T>>;
    friend class Clip<T>;
    constexpr Whole() noexcept: Whole( 0.0, 1.0, 0, 1 ) {}

    constexpr Whole( T lower, T upper, int dim ) noexcept
      : _dim(dim), _delta( (upper - lower) / dim ) {
      _lo_rel = lower / _delta;
    }

    constexpr int dim() const noexcept { return _dim; }
    constexpr T delta() const noexcept { return _delta; }
  };

}

namespace knl :: grid1d {
  template < typename T >
  struct Clip : Expression<Clip<T>> {
  private:
    const Whole<T>& _whole;
    int _anchor; // the cell in the super gridline
    int _extent;

    constexpr T _lower_rel() const noexcept { return _whole._lower_rel() + _anchor; }
  public:
    friend class Expression<Clip<T>>;

    constexpr Clip( const Whole<T>& whole, int anchor, int extent ) noexcept
      : _whole(whole), _anchor(anchor), _extent(extent) {}

    constexpr int dim() const noexcept { return _extent; }
    constexpr T delta() const noexcept { return _whole.delta(); }

    constexpr int anchor() const noexcept { return _anchor; }
    constexpr int extent() const noexcept { return _extent; }

    constexpr void set_anchor( int anchor ) noexcept { _anchor = anchor; }
    constexpr void set_extent( int extent ) noexcept { _extent = extent; }
  };
}

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
