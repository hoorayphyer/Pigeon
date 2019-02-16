#include "particle/pusher.cpp"
#include "traits.hpp"
#include "apt/vec.hpp"
#include "particle/particle.hpp"

using namespace traits;

namespace particle {
#define INSTANTIATE_SPECIES(SP)                                 \
  template apt::Vec<real_t, DPtc>                               \
  update_p< species::SP,                                        \
            apt::Vec<real_t, DPtc>,                             \
            vParticle<real_t,DPtc,ptc_state_t>,                 \
            apt::Vec<real_t, DPtc>,                             \
            real_t >                                            \
  ( PtcExpression<vParticle<real_t, DPtc, ptc_state_t> >& ptc,  \
    const real_t& dt,                                           \
    const apt::VecExpression<apt::Vec<real_t, DPtc>>& E,        \
    const apt::VecExpression<apt::Vec<real_t, DPtc>>& B );      \
                                                                \
  template apt::Vec<real_t, DPtc>                               \
  update_q< species::SP,                                        \
            coordinate_system,                                  \
            apt::Vec<real_t, DPtc>,                             \
            vParticle<real_t,DPtc,ptc_state_t>,                 \
            real_t >                                            \
  ( PtcExpression<vParticle<real_t, DPtc, ptc_state_t> >& ptc,  \
    const real_t& dt )                                          \

  INSTANTIATE_SPECIES(electron);
  INSTANTIATE_SPECIES(positron);
  INSTANTIATE_SPECIES(ion);

}
