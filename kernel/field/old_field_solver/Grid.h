#ifndef _GRID_H_
#define _GRID_H_

#include "Predefs.h"

struct Grid {
  int dims[3]{1, 1, 1};  //!< Dimensions of the grid of each direction
  int guard[3]{};  //!< Number of guard cells at either end of each direction
  int indent[NUM_BOUNDARIES]{};  //!< Indent of the physical domain, required
                                 //!< near the boundaries

  Scalar delta[3]{1.0, 1.0, 1.0};  //!< Grid spacing on each direction (spacing
                                   //!< in coordinate space)
  Scalar lower[3]{};  //!< Lower limit of the grid on each direction

  int dimension = 1;

  ////////////////////////////////////////////////////////////////////////////////
  ///  Reduced dimension in one direction.
  ///
  ///  Reduced dimension means the total size of the grid minus the
  ///  guard cells in both ends. This function is only defined for i >=
  ///  0 and i < DIM.
  ////////////////////////////////////////////////////////////////////////////////
  inline int reducedDim(int i) const { return (dims[i] - 2 * guard[i]); }

  ////////////////////////////////////////////////////////////////////////////////
  ///  Coordinate of a point inside cell n in dimension i.
  ///
  ///  This function applies for both field point and particle
  ///  position. For field, leave pos_in_cell to be default and just
  ///  specify stagger. An unstaggered field component is at the center
  ///  of the cell, while a staggered point is at the end.
  ///
  ///  This calculation is assuming boundary located at the interface
  ///  between guard cells and physical cells. The function is only
  ///  defined for i >= 0 and i < DIM.
  ////////////////////////////////////////////////////////////////////////////////
  inline Scalar pos(int i, int n, int stagger) const {
    // if (i < dimension)
    return (lower[i] + delta[i] * (n - guard[i] + 0.5 + 0.5 * stagger));
    // else
    //  return 0.0;
  }

  inline Scalar pos(int i, int n, StaggerType stagger) const {
    return pos(i, n, static_cast<int>(stagger));
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
  inline Extent extent() const { return Extent{dims[0], dims[1], dims[2]}; }
};

#endif  // ----- #ifndef _GRID_H_  -----
