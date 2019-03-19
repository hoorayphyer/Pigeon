#ifndef _FIELD_SHAPE_INTERPLAY_HPP_
#define _FIELD_SHAPE_INTERPLAY_HPP_

#include "apt/vec.hpp"

namespace field {
  template < typename Field, typename T, typename Vec_q, typename Vec_dq, typename ShapeF >
  void depositWJ ( Field& WJ, T charge,
                   const apt::VecExpression<Vec_q>& q1_abs, // q1 means it's the value after update
                   const apt::VecExpression<Vec_dq>& dq_abs,
                   const ShapeF& shapef );

  template < typename, int, int > struct Field;

  template < typename T, int DField, int DGrid, typename Vec_q, typename ShapeF >
  apt::Vec<T, DField> interpolate ( const Field<T,DField,DGrid>& field,
                                    const apt::VecExpression<Vec_q>& q_abs,
                                    const ShapeF& shapef );

  // the opposite of interpolate. `var` is +=ed to field
  template < typename T, int DField, int DGrid, typename Vec_q, typename ShapeF >
  void deposit ( Field<T,DField,DGrid>& field,
                 apt::Vec<T, DField> var,
                 const apt::VecExpression<Vec_q>& q_abs,
                 const ShapeF& shapef );
}
#endif
