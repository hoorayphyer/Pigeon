#ifndef  _FIELD_COMMUNICATION_HPP_
#define  _FIELD_COMMUNICATION_HPP_

namespace mpi { struct Comm; }

// TODO edge cases: when there is only one node in a dimension, a) what if periodic b) what if not.  Case a) needs some careful treatment because the process may have to send recv within itself, which MPI may not support. Case b) has null neighbors so it shouldn't cause any problem.
namespace field {
  template < typename, int, int > struct Field;

  template < typename T, int DField, int DGrid >
  void sync_guard_cells_from_bulk( Field<T, DField, DGrid>& field, const mpi::CartComm& comm );


  template < typename T, int DField, int DGrid >
  void merge_guard_cells_into_bulk( Field<T, DField, DGrid>& field, const mpi::CartComm& comm );
}

#endif
