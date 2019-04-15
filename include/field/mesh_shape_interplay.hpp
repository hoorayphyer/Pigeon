#ifndef _FIELD_MESH_SHAPE_INTERPLAY_HPP_
#define _FIELD_MESH_SHAPE_INTERPLAY_HPP_

#include "field/field.hpp"
#include "apt/vec.hpp"

namespace field {

  // NOTE q_std refers to the same "standard" as above
  template < typename T, int DField, int DGrid, int Dq, typename ShapeF >
  apt::Vec<T, DField> interpolate ( const Field<T,DField,DGrid>& field,
                                    const apt::array<T,Dq>& q_std,
                                    const ShapeF& shapef );

  // the opposite of interpolate. `var` is +=ed to field // NOTE, this is not same as depositing current
  template < typename T, int DField, int DGrid, int Dq, typename ShapeF >
  void deposit ( Field<T,DField,DGrid>& field,
                 apt::Vec<T, DField> var,
                 const apt::array<T,Dq>& q_std,
                 const ShapeF& shapef );
}
#endif
