#ifndef  _FIELD_COMMUNICATION_HPP_
#define  _FIELD_COMMUNICATION_HPP_

namespace mpi { struct CartComm; }

namespace field {
  template < typename, int, int > struct Field;

  template < typename T, int DField, int DGrid >
  void copy_sync_guard_cells( Field<T, DField, DGrid>& field, const mpi::CartComm& comm );

  template < typename T, int DField, int DGrid >
  void merge_sync_guard_cells( Field<T, DField, DGrid>& field, const mpi::CartComm& comm );
}

#endif
