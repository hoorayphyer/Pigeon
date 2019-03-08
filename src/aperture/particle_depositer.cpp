#include "particle/depositer.cpp"
#include "traits.hpp"

#include "particle/virtual_particle.hpp"
#include "apt/vec.hpp"
#include "apt/grid.hpp"
#include "kernel/shapef.hpp"
using namespace traits;

namespace particle {

  using Ptc = vParticle<real_t, DPtc, ptc_state_t>;
  using Vec = apt::Vec<real_t, DPtc>;
  using Grid = knl::Grid<real_t, DGrid, knl::gridline>;
  using ShapeF = knl::shapef_t<shape>;

  template void
  depositWJ< deposit_j_t, DGrid,
             Ptc,
             Vec,
             Grid,
             ShapeF
             > ( field::Field<deposit_j_t,3,DGrid>& WJ,
                 const PtcExpression<Ptc>& ptc,
                 const apt::VecExpression<Vec>& dq,
                 const Grid& grid,
                 const ShapeF& shapef );
}
