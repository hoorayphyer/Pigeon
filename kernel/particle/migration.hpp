#ifndef  _MIGRATION_HPP_
#define  _MIGRATION_HPP_

#include "apt/vec_expression.hpp"
#include "particle/particle.hpp"
#include "apt/pair.hpp"
#include <vector>

namespace apt {
  constexpr int pow3 ( int i ) noexcept {
    if ( 0 == i ) return 1;
    else return 3 * pow3( i-1 );
  }
}

namespace std { template < class > class optional; }

namespace mpi { struct InterComm; }

namespace particle {
  template < typename Real,
             template < typename > class PtcSpecs,
             int DGrid, int I = DGrid-1 >
  void migrate ( std::vector<Particle< Real, PtcSpecs >>& buffer,
                 const apt::array< apt::pair<std::optional<mpi::InterComm>>, DGrid >& intercomms,
                 unsigned int pairing_shift );
  // NOTE pairing_shift is to help achieve amortized load balancing during intercommunication. This value needs to be synchronized across at least all active processes ( preferrably all processes to avoid future sync. ) Using timestep is a good choice.
}

#endif
