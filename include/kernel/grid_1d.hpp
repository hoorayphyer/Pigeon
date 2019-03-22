#ifndef  _KNL_GRID_1D_HPP_
#define  _KNL_GRID_1D_HPP_

// NOTE convention: the zero of all indices exposed to users start at the first cell in BULK. In other words, guard cells have indices outside [0, dim_of_bulk). Guard is important only when converting from vectorial to linear index, which can be encapsulated in a dedicated function

// guard, indent are controlled by fields directly, where they are collectively called margin cells.

namespace knl :: grid1d {
  template < class E, typename T >
  struct Expression {
  private:
    constexpr T _lower_rel() const noexcept { return static_cast<const E&>(*this)._lower_rel(); }

  public:
    constexpr int dim() const noexcept { return static_cast<const E&>(*this).dim(); }
    constexpr T delta() const noexcept { return static_cast<const E&>(*this).delta(); }

    constexpr T lower() const noexcept { return _lower_rel() * delta(); }
    constexpr T upper() const noexcept { return absc(dim()); }

    // abscissa
    constexpr T absc( int i, T shift_from_lb = 0.0 ) const noexcept {
      return  delta() * ( _lower_rel() + i + shift_from_lb );
    }
  };
}

namespace knl :: grid1d {
  template < typename T > struct Clip;

  template < typename T >
  struct Whole : Expression<Whole<T>, T> {
  private:
    int _dim;
    T _delta;
    T _lo_rel;

    constexpr T _lower_rel() const noexcept { return _lo_rel; }
  public:
    friend class Expression<Whole<T>, T>;
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
  struct Clip : Expression<Clip<T>, T> {
  private:
    const Whole<T>& _whole;
    int _anchor; // the cell in the super gridline
    int _extent;

    constexpr T _lower_rel() const noexcept { return _whole._lower_rel() + _anchor; }
  public:
    friend class Expression<Clip<T>, T>;

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


#endif
