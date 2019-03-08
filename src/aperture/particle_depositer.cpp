#include "particle/depositer.cpp"
#include "traits.hpp"

#include "field/field.hpp"
#include "particle/virtual_particle.hpp"
#include "apt/vec.hpp"
#include "apt/grid.hpp"
#include "kernel/shapef.hpp"
using namespace traits;

namespace particle {

  using Field = field::Field< deposit_j_t, 3, DGrid>;
  using Ptc = vParticle<real_t, DPtc, ptc_state_t>;
  using Vec = apt::Vec<real_t, DPtc>;
  using ShapeF = knl::shapef_t<shape>;

  template void
  depositWJ< Field, Ptc, Vec, ShapeF> ( Field& WJ,
                                              const PtcExpression<Ptc>& ptc,
                                              const apt::VecExpression<Vec>& dq,
                                              const ShapeF& shapef );
}
