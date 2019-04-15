#include "field/mesh_shape_interplay.cpp"
#include "traits.hpp"

#include "field/field.hpp"
#include "apt/vec.hpp"

namespace field {
  using namespace traits;
  constexpr auto Dq = DPtc;

  // NOTE q_std refers to the same "standard" as above
  template apt::Vec<real_t, 3>
  interpolate<real_t, 3, DGrid, Dq, ShapeF> ( const Field<real_t,3,DGrid>& field,
                                          const apt::array<real_t,Dq>& q_std,
                                          const ShapeF& shapef );

  template apt::Vec<real_t, 1>
  interpolate<real_t, 1, DGrid, Dq, ShapeF> ( const Field<real_t, 1, DGrid>& field,
                                          const apt::array<real_t,Dq>& q_std,
                                          const ShapeF& shapef );

  template void
  deposit<real_t, 3, DGrid, Dq, ShapeF> ( Field<real_t, 3, DGrid>& field,
                                      apt::Vec<real_t, 3> var,
                                      const apt::array<real_t,Dq>& q_std,
                                      const ShapeF& shapef );

  template void
  deposit<real_t, 1, DGrid, Dq, ShapeF> ( Field<real_t, 1, DGrid>& field,
                                      apt::Vec<real_t, 1> var,
                                      const apt::array<real_t,Dq>& q_std,
                                      const ShapeF& shapef );
}
