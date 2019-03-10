#ifndef _FIELD_SHAPE_INTERPLAY_HPP_
#define _FIELD_SHAPE_INTERPLAY_HPP_

#include "apt/vec_expression.hpp"

// TODO check missing 4pi's maybe
namespace field {
  // NOTE TWJ here may use long double
  template < typename Field,
             typename Vec_q,
             typename Vec_dq,
             typename ShapeF >
  void depositWJ ( Field& WJ,
                   const apt::VecExpression<Vec_q>& q1_abs, // q1 means it's the value after update
                   const apt::VecExpression<Vec_dq>& dq_abs,
                   const ShapeF& shapef );
}

#endif
