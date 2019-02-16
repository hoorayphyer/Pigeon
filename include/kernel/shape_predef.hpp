#ifndef  _SHAPE_PREDEF_HPP_
#define  _SHAPE_PREDEF_HPP_

// value here denotes length of support
namespace knl {
  enum class shape : int {
    Nearest_Grid_Point = 1,
    Cloud_In_Cell = 2,
    Triangular_Cloud = 3,
    Piecewise_Cubic_Spline = 4
  };
}

#endif
