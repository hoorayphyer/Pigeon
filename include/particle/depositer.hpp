#ifndef _DEPOSITER_HPP_
#define _DEPOSITER_HPP_

// TODO check missing 4pi's maybe
namespace particle {
  // NOTE TWJ here may use long double
  template < typename Field, typename Ptc, typename Vec, typename Grid, typename ShapeF >
  void depositWJ ( Field& WJ, const Ptc& ptc, const Vec& dq, const Grid& grid, const ShapeF& shapef );
}

#endif
