#include "field/communication.cpp"
#include "traits.hpp"

using namespace traits;

namespace field {
  template
  void sync_guard_cells_from_bulk<real_t, 1, DGrid>( Field<real_t, 1, DGrid>& field, const mpi::CartComm& comm );

  template
  void sync_guard_cells_from_bulk<real_t, 3, DGrid>( Field<real_t, 3, DGrid>& field, const mpi::CartComm& comm );

  template
  void merge_guard_cells_into_bulk<real_dj_t, 3, DGrid>( Field<real_dj_t, 3, DGrid>& field, const mpi::CartComm& comm );

}
