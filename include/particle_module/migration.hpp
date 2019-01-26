#ifndef  _MIGRATION_HPP_
#define  _MIGRATION_HPP_

#include "core/particle_vector.hpp"

// TODO may use intercommunicator
namespace mpi {
  struct Comm;
}

namespace particle {
  template < typename Tvt, std::size_t DGrid, typename Trl = apt::remove_cvref_t<Tvt> >
  bool is_migrate( const Vec<Tvt, DGrid>& q, const std::array< std::array<Trl, 2>, DGrid>& bounds ) noexcept;

  template < std::size_t DGrid, typename T, std::size_t DPtc, std::size_t I = DGrid - 1 >
  void migrate( particle::vector<T,DPtc>& buffer,
                const std::array< std::array<int,2>, DGrid >& neighbors,
                const std::array< std::array<T,2>, DGrid>& bounds,
                const mpi::Comm& comm );

}

#endif
