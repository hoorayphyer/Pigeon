#include "field/field_shape_interplay.cpp"
#include "traits.hpp"

#include "field/field.hpp"
#include "apt/vec.hpp"
#include "apt/grid.hpp"
#include "kernel/shapef.hpp"
using namespace traits;

namespace field {

  using Field = Field< deposit_j_t, 3, DGrid>;
  using Vec_q = apt::vVec<real_t, DPtc>;
  using Vec_dq = apt::Vec<real_t, DPtc>;
  using ShapeF = knl::shapef_t<shape>;

  template void
  depositWJ< Field, Vec_q, Vec_dq, ShapeF> ( Field& WJ,
                                             const apt::VecExpression<Vec_q>& q1_abs,
                                             const apt::VecExpression<Vec_dq>& dq_abs,
                                             const ShapeF& shapef );
}
