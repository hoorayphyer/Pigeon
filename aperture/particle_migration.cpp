#include "particle/migration.cpp"
#include "traits.hpp"

#include "particle/virtual_particle.hpp"

using namespace traits;

namespace particle {
  template void
  migrate < real_t, DPtc, ptc_state_t, DGrid >
  ( std::vector<cParticle<real_t,DPtc,ptc_state_t>>& buffer,
    const apt::array< apt::pair<std::optional<mpi::InterComm>>, DGrid >& intercomms,
    unsigned int pairing_shift );
}
