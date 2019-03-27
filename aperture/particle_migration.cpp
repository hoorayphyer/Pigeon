#include "particle/migration.cpp"
#include "traits.hpp"

#include "particle/virtual_particle.hpp"

using namespace traits;

namespace particle {
  template bool
  is_migrate< apt::vVec<real_t,DPtc>, DGrid, real_t >
  ( const apt::VecExpression<apt::vVec<real_t,DPtc>, real_t>& q,
    const apt::array< apt::pair<real_t>, DGrid>& borders ) noexcept;

  template void
  migrate < real_t, DPtc, ptc_state_t, DGrid >
  ( std::vector<cParticle<real_t,DPtc,ptc_state_t>>& buffer,
    const apt::array< apt::pair<std::optional<mpi::InterComm>>, DGrid >& intercomms,
    const apt::array< apt::pair<real_t>, DGrid>& borders,
    unsigned int pairing_shift );
}
