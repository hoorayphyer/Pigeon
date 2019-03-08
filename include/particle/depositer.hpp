#ifndef _DEPOSITER_HPP_
#define _DEPOSITER_HPP_

#include "particle/particle_expression.hpp"
#include "apt/vec_expression.hpp"

// TODO check missing 4pi's maybe
namespace particle {
  // NOTE TWJ here may use long double
  template < typename Field,
             typename Ptc,
             typename Vec,
             typename ShapeF >
  void depositWJ ( Field& WJ,
                   const PtcExpression<Ptc>& ptc,
                   const apt::VecExpression<Vec>& dq,
                   const ShapeF& shapef );
}

#endif
