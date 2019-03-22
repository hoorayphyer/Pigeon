#include "particle/pusher.cpp"
#include "traits.hpp"
#include "apt/vec.hpp"
#include "particle/virtual_particle.hpp"

using namespace traits;

namespace particle {
  template apt::Vec<real_t, DPtc>
  update_p< apt::Vec<real_t, DPtc>,
            vParticle<real_t,DPtc,ptc_state_t>,
            apt::Vec<real_t, DPtc>,
            real_t >
  ( PtcExpression<vParticle<real_t, DPtc, ptc_state_t> >& ptc,
    real_t dt, unsigned int mass_x,
    const apt::VecExpression<apt::Vec<real_t, DPtc>>& E,
    const apt::VecExpression<apt::Vec<real_t, DPtc>>& B );

  template apt::Vec<real_t, DPtc>
  update_q< coordinate_system,
            apt::Vec<real_t, DPtc>,
            vParticle<real_t,DPtc,ptc_state_t>,
            real_t >
  ( PtcExpression<vParticle<real_t, DPtc, ptc_state_t> >& ptc,
    real_t dt, bool is_massive );
}
