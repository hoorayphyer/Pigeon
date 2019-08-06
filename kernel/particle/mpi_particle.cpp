#include "particle/mpi_particle.hpp"
#include "mpipp/mpi_p2p_impl.hpp"
#include "mpipp/mpi_collective_impl.hpp"
#include "pic.hpp"

namespace mpi {
  using namespace pic;

  using Ptc_t = particle::Particle<real_t, particle::Specs>;
  template struct Datatype<Ptc_t>;
  INSTANTIATE_MPI_P2P(Ptc_t);
  INSTANTIATE_MPI_COLLECTIVE(Ptc_t);
}
