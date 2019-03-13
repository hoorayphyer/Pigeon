#ifndef _GRID_H_
#define _GRID_H_

#include "Predefs.h"

// template <int DIM>
struct Grid {
  int dims[3];   //!< Dimensions of the grid of each direction
  int guard[3];  //!< Number of guard cells at either end of each direction
  int indent[NUM_BOUNDARIES]; //!< Indent of the physical domain, required near the boundaries

  Scalar delta[3];  //!< Grid spacing on each direction (spacing in coordinate
                    //!< space)
  Scalar lower[3];  //!< Lower limit of the grid on each direction

  int dimension;

  Grid() {
    // Initialize all quantities to zero, and dimensions to 1
    for (int i = 0; i < 3; i++) {
      dims[i] = 1;
      guard[i] = 0;
      delta[i] = 1.0;
      lower[i] = 0.0;
    }
    dimension = 1;
  }

  ////////////////////////////////////////////////////////////////////////////////
  ///  Constructor which only initialize dimensions.
  ////////////////////////////////////////////////////////////////////////////////
  Grid(int N1, int N2 = 1, int N3 = 1) {
    dims[0] = N1;
    dims[1] = N2;
    dims[2] = N3;

    // Initialize other quantities to zero
    for (int i = 0; i < 3; i++) {
      guard[i] = 0;
      delta[i] = 1.0;
      lower[i] = 0.0;
    }

    if (dims[1] <= 1 && dims[2] <= 1)
      dimension = 1;
    else if (dims[2] <= 1)
      dimension = 2;
    else
      dimension = 3;
  }

  ////////////////////////////////////////////////////////////////////////////////
  ///  Reduced dimension in one direction.
  ///
  ///  Reduced dimension means the total size of the grid minus the
  ///  guard cells in both ends. This function is only defined for i >=
  ///  0 and i < DIM.
  ////////////////////////////////////////////////////////////////////////////////
  // template <int i,
  //           typename = typename std::enable_if<(i >= 0 && i < DIM)>::type>
  inline int reducedDim(int i) const { return (dims[i] - 2 * guard[i]); }

  ////////////////////////////////////////////////////////////////////////////////
  ///  Coordinate of a point inside cell n in dimension i.
  ///
  ///  This function applies for both field point and particle
  ///  position. For field, leave pos_in_cell to be default and just
  ///  specify stagger. An unstaggered field component is at the center
  ///  of the cell, while a staggered point is at the end. For
  ///  particles, specify pos_in_cell as the relative position of the
  ///  particle in the cell and put stagger equal to zero.
  ///
  ///  This calculation is assuming boundary located at the interface
  ///  between guard cells and physical cells. The function is only
  ///  defined for i >= 0 and i < DIM.
  ////////////////////////////////////////////////////////////////////////////////
  // template <int i,
  //           typename = typename std::enable_if<(i >= 0 && i < DIM)>::type>
  inline Scalar pos(int i, int n, int stagger, Scalar pos_in_cell = 0.5) const {
   // if (i < dimension)
      return (lower[i] + delta[i] * (n - guard[i] + pos_in_cell + 0.5 * stagger));
   // else
    //  return 0.0;
  }

  inline Scalar pos(int i, int n, StaggerType stagger,
             Scalar pos_in_cell = 0.5) const {
    //if (i < dimension)
      return (lower[i] +
              delta[i] *
              (n - guard[i] + pos_in_cell + 0.5 * static_cast<int>(stagger)));
    //else
     // return 0.0;
  }

  inline Vec3<Scalar> pos_3d(int idx, const Vec3<int>& stagger) const {
    Vec3<Scalar> result;
    result[0] = pos(0, idx % dims[0], stagger.x);
    result[1] = pos(1, (idx / dims[0]) % dims[1], stagger.y);
    result[2] = pos(2, idx / (dims[0] * dims[1]), stagger.z);
    return result;
  }

  // given a scalar position q, find the cell n in that direction, such that q =
  // lower + delta * ( n - guard + x ), where 0 <= x < 1.
  inline int pos_to_cell( Scalar pos, int dir ) const {
    return static_cast<int>( ( pos - lower[dir] ) / delta[dir] );
  }

  ////////////////////////////////////////////////////////////////////////////////
  ///  Full coordinate of a point inside the grid.
  ///
  ///  This function applies for only particle position. Pos_rel has
  ///  components greater than zero and less than the corresponding
  ///  delta.
  ///
  ///  This calculation is assuming boundary located at the interface
  ///  between guard cells and physical cells. The function is only
  ///  defined for i >= 0 and i < DIM.
  ////////////////////////////////////////////////////////////////////////////////
  template <typename T>
  inline Vec3<Scalar> pos_particle( int cell_linear, const Vec3<T>& pos_rel ) const {
      Vec3<Scalar> pos_full( pos(0, getC1(cell_linear), 0, pos_rel.x),
                             pos(1, getC2(cell_linear), 0, pos_rel.y),
                             pos(2, getC3(cell_linear), 0, pos_rel.z)
                             );//Note deltas cannot be zero, even that dimension is absent.
      return pos_full;
  }

  ////////////////////////////////////////////////////////////////////////////////
  ///  Upper boundary position in direction i
  ////////////////////////////////////////////////////////////////////////////////
  inline Scalar upper( int i ) const {
    return pos(i, dims[i] - guard[i] - 1, 1);
  }

  ////////////////////////////////////////////////////////////////////////////////
  ///  Index of the point if the grid is stratified into 1 direction.
  ////////////////////////////////////////////////////////////////////////////////
  inline int getIdx(int c1, int c2 = 0, int c3 = 0) const {
    return c1 + c2 * dims[0] + c3 * dims[0] * dims[1];
  }

  ////////////////////////////////////////////////////////////////////////////////
  ///  Test if a point is inside the grid.
  ////////////////////////////////////////////////////////////////////////////////
  inline bool isInGrid(int c1, int c2 = 0, int c3 = 0) const {
    return (c1 >= 0 && c1 < dims[0]) && (c2 >= 0 && c2 < dims[1]) &&
           (c3 >= 0 && c3 < dims[2]);
  }

  ////////////////////////////////////////////////////////////////////////////////
  ///  Test if a point is inside the bulk of the grid, not in guard cells.
  ////////////////////////////////////////////////////////////////////////////////
  inline bool isInBulk(int c1, int c2 = 0, int c3 = 0) const {
    return (c1 >= guard[0] && c1 < dims[0] - guard[0])
        && (c2 >= guard[1] && c2 < dims[1] - guard[1])
        && (c3 >= guard[2] && c3 < dims[2] - guard[2]);
  }

  ////////////////////////////////////////////////////////////////////////////////
  ///  Test if a point is inside the bulk of the grid, not in guard cells.
  ////////////////////////////////////////////////////////////////////////////////
  inline bool isInBulk(const Index& idx) const {
    return (idx.x >= guard[0] && idx.x < dims[0] - guard[0])
        && (idx.y >= guard[1] && idx.y < dims[1] - guard[1])
        && (idx.z >= guard[2] && idx.z < dims[2] - guard[2]);
  }

  ////////////////////////////////////////////////////////////////////////////////
  ///  Get the size of the grid (product of all dimensions).
  ////////////////////////////////////////////////////////////////////////////////
  inline int size() const {
    int tmp = 1;
    // #pragma unroll
    for (int i = 0; i < 3; i++) {
      tmp *= dims[i];
    }
    return tmp;
  }

  ////////////////////////////////////////////////////////////////////////////////
  ///  Get the extent of the grid. Used for interfacing with multiarrays.
  ////////////////////////////////////////////////////////////////////////////////
  inline Extent extent() const {
    return Extent{dims[0], dims[1], dims[2]};
    //    return tmp;
  }

  inline Extent extent_less() const {
    return Extent{dims[0] - 2 * guard[0], dims[1] - 2 * guard[1],
                  dims[2] - 2 * guard[2]};
  }

  inline int getC1(int idx) const { return idx % dims[0]; }
  inline int getC2(int idx) const { return (idx / dims[0]) % dims[1]; }
  inline int getC3(int idx) const { return idx / (dims[0] * dims[1]); }
  inline Vec3<int> getCell(int idx) const {
    return Vec3<int>( getC1(idx), getC2(idx), getC3(idx) );
  }

};

#endif  // ----- #ifndef _GRID_H_  -----
