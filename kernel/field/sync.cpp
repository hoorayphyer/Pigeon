#include "field/sync_impl.hpp"
#include "pic.hpp"

using namespace pic;

// TODO If real_t = real_j_t or real_export_t, there will be duplicate instantiation. For now we use just use explicit types
namespace field {
  template
  void copy_sync_guard_cells( Field<float, 1, DGrid>& , const mpi::CartComm&,
                              const apt::array<apt::Range,DGrid>& domain_range,
                              const apt::array<apt::Range,DGrid>& sub_range );
  template
  void copy_sync_guard_cells( Field<float, 3, DGrid>& , const mpi::CartComm&,
                              const apt::array<apt::Range,DGrid>& domain_range,
                              const apt::array<apt::Range,DGrid>& sub_range );
  template
  void merge_sync_guard_cells( Field<float, 3, DGrid>& , const mpi::CartComm&,
                               const apt::array<apt::Range,DGrid>& domain_range,
                               const apt::array<apt::Range,DGrid>& sub_range );
  template
  void merge_sync_guard_cells( Field<float, 1, DGrid>& , const mpi::CartComm&,
                               const apt::array<apt::Range,DGrid>& domain_range,
                               const apt::array<apt::Range,DGrid>& sub_range );

  template
  void copy_sync_guard_cells( Field<double, 1, DGrid>& , const mpi::CartComm&,
                              const apt::array<apt::Range,DGrid>& domain_range,
                              const apt::array<apt::Range,DGrid>& sub_range );
  template
  void copy_sync_guard_cells( Field<double, 3, DGrid>& , const mpi::CartComm&,
                              const apt::array<apt::Range,DGrid>& domain_range,
                              const apt::array<apt::Range,DGrid>& sub_range );
  template
  void merge_sync_guard_cells( Field<double, 3, DGrid>& , const mpi::CartComm&,
                               const apt::array<apt::Range,DGrid>& domain_range,
                               const apt::array<apt::Range,DGrid>& sub_range );
  template
  void merge_sync_guard_cells( Field<double, 1, DGrid>& , const mpi::CartComm&,
                               const apt::array<apt::Range,DGrid>& domain_range,
                               const apt::array<apt::Range,DGrid>& sub_range );
}
