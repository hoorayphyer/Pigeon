#include "particle/annihilation_impl.hpp"
#include "pic.hpp"

using namespace pic;

namespace particle {
  template void annihilate<DGrid, real_t, Specs, ShapeF, real_j_t>
  ( array<real_t,Specs>& el, array<real_t,Specs>& po,
    field::Field<real_j_t,3,DGrid>& J,
    real_t charge_el, real_t charge_po,
    const apt::Grid< real_t, DGrid >& grid,
    const mpi::Comm& intra,
    real_t dt, const ShapeF&,
    real_t(*)(real_t, real_t) );
}
