#ifndef  _MIGRATION_HPP_
#define  _MIGRATION_HPP_

#include "apt/vec_expression.hpp"
#include "particle/c_particle.hpp"
#include "apt/pair.hpp"
#include <vector>

namespace std { template < class > class optional; }

namespace mpi { struct InterComm; }

namespace particle {
  template < typename Vec, int DGrid, typename T >
  bool is_migrate( const apt::VecExpression<Vec,T>& q,
                   const apt::array< apt::pair<T>, DGrid>& borders ) noexcept;

  template < typename T, int DPtc, typename state_t, int DGrid >
  void migrate ( std::vector<cParticle<T,DPtc,state_t>>& buffer,
                 const apt::array< apt::pair<std::optional<mpi::InterComm>>, DGrid >& intercomms,
                 const apt::array< apt::pair<T>, DGrid>& borders,
                 unsigned int pairing_shift );
  // NOTE pairing_shift is to help achieve amortized load balancing during intercommunication. This value needs to be synchronized across at least all active processes ( preferrably all processes to avoid future sync. ) Using timestep is a good choice.
}

#endif
