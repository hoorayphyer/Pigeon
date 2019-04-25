#include "field/mesh_shape_interplay.cpp"
#include "pic.hpp"

namespace field {
  using namespace pic;
  constexpr auto Dq = particle::Specs<real_t>::Dim;

  // NOTE q_std refers to the same "standard" as above
  template apt::Vec<real_t, 3>
  interpolate ( const Field<real_t,3,DGrid>& field,
                const apt::array<real_t,Dq>& q_std,
                const ShapeF& shapef );

  template apt::Vec<real_t, 1>
  interpolate ( const Field<real_t, 1, DGrid>& field,
                const apt::array<real_t,Dq>& q_std,
                const ShapeF& shapef );

  template apt::Vec<real_t, 3>
  interpolate ( const Field<real_t,3,DGrid>& field,
                const apt::array<real_t,DGrid>& q_std, // NOTE it is DGrid instead of Dq. This is version is used in dataexport
                const ShapeF& shapef );

  template apt::Vec<real_t, 1>
  interpolate ( const Field<real_t, 1, DGrid>& field,
                const apt::array<real_t,DGrid>& q_std, // NOTE it is DGrid instead of Dq.
                const ShapeF& shapef );

  template void
  deposit ( Field<real_t, 3, DGrid>& field,
            apt::array<real_t, 3> var,
            const apt::array<real_t,Dq>& q_std,
            const ShapeF& shapef );

  template void
  deposit ( Field<real_t, 1, DGrid>& field,
            apt::array<real_t, 1> var,
            const apt::array<real_t,Dq>& q_std,
            const ShapeF& shapef );
}
