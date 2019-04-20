#ifndef _FIELD_CURRENT_DEPOSITION_HPP_
#define _FIELD_CURRENT_DEPOSITION_HPP_

#include "field/field.hpp"
#include "apt/array.hpp"

namespace field {
  // NOTE The standard grid of the bulk of the mesh is a rescaled and shifted version of the actual grid such that the spacing = 1 and the first cell has index = 0. This way, there is no need of grid information. "Standard" is borrowed from "standard normal distribution".
  template < typename RealJ, int DField, int DGrid, typename ShapeF, typename U >
  void deposit ( Field<RealJ,DField,DGrid>& J,
                 U charge_over_dt,
                 const ShapeF& shapef,
                 const apt::array<U,DField>& q0_std, // NOTE its DField, not DGrid
                 const apt::array<U,DField>& q1_std );

  template < typename RealJ, int DField, int DGrid >
  void integrate( Field<RealJ,DField,DGrid>& J );


}

#endif
