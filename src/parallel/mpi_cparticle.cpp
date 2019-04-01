#ifndef  _MPI_DATATYPE_CPARTICLE_HPP_
#define  _MPI_DATATYPE_CPARTICLE_HPP_
#include "parallel/mpi_datatype.hpp"

#include "particle/c_particle.hpp"
#include "../../aperture/traits.hpp"
// TODOL use typelist.
// See also mpi++.cpp for type commit
// See also main.cpp for type commit
namespace mpi {
  using namespace traits;
  using cPtc_t = particle::cParticle<real_t,DPtc,ptc_state_t>;
}
#endif
