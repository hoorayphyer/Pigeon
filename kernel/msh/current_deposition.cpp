#include "msh/current_deposition_impl.hpp"
#include "pic.hpp"

namespace msh {
using namespace pic;

template void deposit<real_j_t, 3, DGrid, ShapeF, real_t>(
    field::Field<real_j_t, 3, DGrid>&, real_t, const ShapeF&,
    const apt::array<real_t, 3>&, const apt::array<real_t, 3>&);

// template
// void integrate< real_j_t, 3, DGrid > ( field::Field<real_j_t, 3, DGrid>& );
}  // namespace msh
