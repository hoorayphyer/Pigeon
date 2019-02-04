#ifndef  _MIGRATION_HPP_
#define  _MIGRATION_HPP_

#include <array>

namespace mpi { struct Comm; }


// TODO may use intercommunicator

namespace particle {
  template < typename Vec, int DGrid, typename T >
  bool is_migrate( const Vec& q, const std::array< std::array<T, 2>, DGrid>& borders ) noexcept;

  template < typename PtcArray, int DGrid, typename T >
  void migrate ( PtcArray& buffer, const std::array< std::array<int,2>, DGrid >& neighbors,
                 const std::array< std::array<T,2>, DGrid>& borders, const mpi::Comm& comm );
}

#endif
