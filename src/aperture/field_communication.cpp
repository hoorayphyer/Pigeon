#include "field/communication.cpp"
#include "traits.hpp"

using namespace traits;

namespace field {
  template
  void sync_guard_cells<real_t, 1, DGrid>( Field<real_t, 1, DGrid>& field, const mpi::Comm& comm );

  template
  void sync_guard_cells<real_t, 3, DGrid>( Field<real_t, 3, DGrid>& field, const mpi::Comm& comm );

}
