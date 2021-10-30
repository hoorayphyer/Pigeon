#ifndef _APT_GRID_HPP_
#define _APT_GRID_HPP_

#include <cmath>

#include "apt/array.hpp"

// NOTE convention: the zero of all indices exposed to users start at the first
// cell in BULK. In other words, guard cells have indices outside [0,
// dim_of_bulk). Guard is important only when converting from vectorial to
// linear index, which can be encapsulated in a dedicated function

// guard, indent are controlled by fields directly, where they are collectively
// called margin cells.

// namespace knl {
//   template < typename T >
//   struct SuperGrid1D {
//   private:
//     int _dim;
//     T _delta;
//     T _lower;

//   public:
//     constexpr SuperGrid1D() noexcept: SuperGrid1D( 0.0, 1.0, 1 ) {}

//     constexpr SuperGrid1D( T lower, T upper, int dim ) noexcept
//       : _dim(dim),
//         _delta( (upper - lower) / dim ),
//         _lower(lower) {}

//     constexpr int dim() const noexcept { return _dim; }
//     constexpr T delta() const noexcept { return _delta; }

//     constexpr T lower() const noexcept { return _lower; }
//     constexpr T upper() const noexcept { return absc(dim()); }

//     // abscissa
//     constexpr T absc( int i, T shift_from_lb = 0.0 ) const noexcept {
//       return  _lower + delta() * ( i + shift_from_lb );
//     }
//   };

// }

// namespace knl {
//   template < typename T >
//   struct SubGrid1D {
//   private:
//     const SuperGrid1D<T>& _supergrid;
//     int _anchor; // the cell in the super gridline
//     int _dim;

//   public:
//     constexpr SubGrid1D( const SuperGrid1D<T>& supergrid, int anchor, int dim
//     ) noexcept
//       : _supergrid(supergrid), _anchor(anchor), _dim(dim) {}

//     constexpr int dim() const noexcept { return _dim; }
//     constexpr T delta() const noexcept { return _supergrid.delta(); }

//     constexpr T lower() const noexcept { return _supergrid.absc(_anchor); }
//     constexpr T upper() const noexcept { return _supergrid.absc(_anchor +
//     _dim); }
//   };
// }

namespace apt {
template <typename T>
struct Grid1D {
 private:
  int _dim;
  T _delta;
  T _lower;

 public:
  constexpr Grid1D() noexcept : Grid1D(0.0, 1.0, 1) {}

  constexpr Grid1D(T lower, T upper, int dim) noexcept
      : _dim(dim), _delta((upper - lower) / dim), _lower(lower) {}

  constexpr int dim() const noexcept { return _dim; }
  constexpr T delta() const noexcept { return _delta; }

  constexpr T lower() const noexcept { return _lower; }
  constexpr T upper() const noexcept { return absc(dim()); }

  // abscissa
  constexpr T absc(int i, T shift_from_lb = 0.0) const noexcept {
    return _lower + delta() * (i + shift_from_lb);
  }

  // reverse abscissa, find index from abscissa
  constexpr int csba(T x) const noexcept {
    // TODO optimize it?
    return (x < _lower)
               ? -static_cast<int>(std::ceil((_lower - x) / delta()) + 0.5)
               : static_cast<int>((x - _lower) / delta());
  }

  constexpr void clip(int i_start, int extent) noexcept {
    _lower = absc(i_start);
    _dim = extent;
  }

  constexpr Grid1D divide(int num_pieces, int ith_piece) const noexcept {
    Grid1D res = *this;
    int dim = _dim / num_pieces;
    res.clip(ith_piece * dim, dim);
    return res;
  }
};

}  // namespace apt

namespace apt {
template <typename T, int DGrid>
using Grid = array<Grid1D<T>, DGrid>;

template <typename T, int DGrid>
constexpr array<int, DGrid> dims(const Grid<T, DGrid>& grid) noexcept {
  array<int, DGrid> ext;
  for (int i = 0; i < DGrid; ++i) ext[i] = grid[i].dim();
  return ext;
}

template <typename T, int DGrid>
constexpr array<T, DGrid> abscs(const Grid<T, DGrid>& grid,
                                const array<int, DGrid>& I,
                                const array<T, DGrid>& shift = {}) noexcept {
  array<T, DGrid> res;
  for (int i = 0; i < DGrid; ++i) res[i] = grid[i].absc(I[i], shift[i]);
  return res;
}

// promote the dimension
template <int D, typename T, int DGrid>
constexpr array<T, D> abscs(const Grid<T, DGrid>& grid,
                            const array<int, DGrid>& I,
                            const array<T, DGrid>& shift = {}) noexcept {
  static_assert(D >= DGrid);
  array<T, D> res{};
  for (int i = 0; i < DGrid; ++i) res[i] = grid[i].absc(I[i], shift[i]);
  return res;
}

template <typename T, int DGrid>
constexpr T dV(const Grid<T, DGrid>& grid) noexcept {
  T res = 1.0;
  for (const auto& g : grid) res *= g.delta();
  return res;
}
}  // namespace apt

#endif
