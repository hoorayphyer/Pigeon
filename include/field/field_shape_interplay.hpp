#ifndef _FIELD_SHAPE_INTERPLAY_HPP_
#define _FIELD_SHAPE_INTERPLAY_HPP_

#include "apt/vec.hpp"

namespace field {
  template < typename Field, typename T, typename Vec_q, typename Vec_dq, typename ShapeF >
  void deposit_dJ ( Field& dJ, T charge,
                    const apt::VecExpression<Vec_q>& q0_abs, // q1 = q0 + dq
                    const apt::VecExpression<Vec_dq>& dq_abs,
                    const ShapeF& shapef );

  template < typename, int, int > struct Field;

  template < typename T, int DField, int DGrid, typename LocType, typename ShapeF >
  apt::Vec<T, DField> interpolate ( const Field<T,DField,DGrid>& field,
                                    const LocType& q_abs,
                                    const ShapeF& shapef );

  // the opposite of interpolate. `var` is +=ed to field
  template < typename T, int DField, int DGrid, typename LocType, typename ShapeF >
  void deposit ( Field<T,DField,DGrid>& field,
                 apt::Vec<T, DField> var,
                 const LocType& q_abs,
                 const ShapeF& shapef );
}
#endif
