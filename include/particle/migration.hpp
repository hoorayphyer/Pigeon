#ifndef  _MIGRATION_HPP_
#define  _MIGRATION_HPP_

#include "apt/vec_expression.hpp"
#include <vector>

namespace std {
  template < typename, std::size_t > struct array;
  template < class > class optional;
}

namespace mpi { struct InterComm; }
namespace particle {
  template <typename,int,typename> struct cParticle;
}

namespace particle {
  template < typename Vec, int DGrid, typename T >
  bool is_migrate( const apt::VecExpression<Vec>& q,
                   const std::array< std::array<T, 2>, DGrid>& borders ) noexcept;

  template < typename T, int DPtc, typename state_t, int DGrid >
  void migrate ( std::vector<cParticle<T,DPtc,state_t>>& buffer,
                 const std::array< std::array<std::optional<mpi::InterComm>,2>, DGrid >& intercomms,
                 const std::array< std::array<T,2>, DGrid>& borders, unsigned int pairing_shift );
  // NOTE pairing_shift is to help achieve amortized load balancing during intercommunication. This value needs to be synchronized across at least all active processes ( preferrably all processes to avoid future sync. ) Using timestep is a good choice.
}

#endif
