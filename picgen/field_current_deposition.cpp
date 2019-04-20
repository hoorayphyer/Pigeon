#include "field/current_deposition.cpp"
#include "pic.hpp"

namespace field {
  using namespace pic;

  template
  void deposit<real_j_t, 3, DGrid, ShapeF, real_t> ( Field<real_j_t,3,DGrid>&, real_t, const ShapeF&, const apt::array<real_t,3>&, const apt::array<real_t,3>& );

  template
  void integrate< real_j_t, 3, DGrid > ( Field<real_j_t, 3, DGrid>& );
}
