#ifndef _DEPOSITER_HPP_
#define _DEPOSITER_HPP_

#include "particle/particle_expression.hpp"
#include "apt/vec_expression.hpp"

namespace knl { template < int, typename > struct Grid; }
namespace field { template < typename, int, int > struct Field; }

// TODO check missing 4pi's maybe
namespace particle {
  // NOTE TWJ here may use long double
  template < typename T_deposit_j, int DGrid,
             typename Ptc,
             typename Vec,
             typename T,
             typename ShapeF >
  void depositWJ ( field::Field<T_deposit_j,3,DGrid>& WJ,
                   const PtcExpression<Ptc>& ptc,
                   const apt::VecExpression<Vec>& dq,
                   const knl::Grid<DGrid,T>& grid,
                   const ShapeF& shapef );
}

#endif
