#include "field/communication.cpp"
#include "traits.hpp"

using namespace traits;

namespace field {
  template
  void sync_guard_cells_from_bulk( Field<real_t, 1, DGrid>& field, const mpi::CartComm& comm );

  template
  void sync_guard_cells_from_bulk( Field<real_t, 3, DGrid>& field, const mpi::CartComm& comm );

  template
  void merge_guard_cells_into_bulk( Field<real_j_t, 3, DGrid>& field, const mpi::CartComm& comm );

}
