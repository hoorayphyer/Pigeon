#include "particle/migration.cpp"
#include "traits.hpp"

#include "particle/virtual_particle.hpp"

using namespace traits;

namespace particle {
  template bool
  is_migrate< apt::vVec<real_t,DPtc>,
              DGrid,
              real_t > ( const apt::VecExpression<apt::vVec<real_t,DPtc>>& q,
                         const std::array< std::array<real_t, 2>, DGrid>& borders ) noexcept;

  template void
  migrate < real_t, DPtc, ptc_state_t, DGrid >
  ( std::vector<cParticle<real_t,DPtc,ptc_state_t>>& buffer,
    const std::array< std::array<std::optional<mpi::InterComm>,2>, DGrid >& intercomms,
    const std::array< std::array<real_t,2>, DGrid>& borders,
    unsigned int pairing_shift );
}
