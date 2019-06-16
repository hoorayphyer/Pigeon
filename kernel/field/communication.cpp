#include "communication_impl.hpp"
#include "pic.hpp"

using namespace pic;

// TODO what if real_t = real_j_t or real_export_t
namespace field {
  template
  void sync_guard_cells_from_bulk( Field<real_t, 1, DGrid>& field, const mpi::CartComm& comm );

  template
  void sync_guard_cells_from_bulk( Field<real_t, 3, DGrid>& field, const mpi::CartComm& comm );

  template
  void sync_guard_cells_from_bulk( Field<real_j_t, 3, DGrid>& field, const mpi::CartComm& comm );

  template
  void merge_guard_cells_into_bulk( Field<real_j_t, 3, DGrid>& field, const mpi::CartComm& comm );

}

namespace field {
  template
  void sync_guard_cells_from_bulk( Field<real_export_t, 1, DGrid>& field, const mpi::CartComm& comm );
}
