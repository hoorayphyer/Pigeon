#ifndef  _MPI_DATATYPE_PARTICLE_HPP_
#define  _MPI_DATATYPE_PARTICLE_HPP_
#include "mpipp/mpi_datatype.hpp"

#include "particle/particle.hpp"
#include "pic.hpp"
// TODOL use typelist.
// See also mpi++.cpp for type commit
namespace mpi {
  using namespace pic;
  using Ptc_t = particle::Particle<real_t, particle::Specs>;
}
#endif
