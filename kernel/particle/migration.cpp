#include "migration_impl.hpp"
#include "pic.hpp"

using namespace pic;

namespace particle {
  template void migrate ( std::vector<cParticle< real_t, Specs >>& buffer,
                          const apt::array< apt::pair<std::optional<mpi::InterComm>>, DGrid >& intercomms,
                          unsigned int pairing_shift );
}
