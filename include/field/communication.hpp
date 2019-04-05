#ifndef  _FIELD_COMMUNICATION_HPP_
#define  _FIELD_COMMUNICATION_HPP_

namespace mpi { struct CartComm; }

namespace field {
  template < typename, int, int > struct Field;

  template < typename T, int DField, int DGrid >
  void sync_guard_cells_from_bulk( Field<T, DField, DGrid>& field, const mpi::CartComm& comm );


  template < typename T, int DField, int DGrid >
  void merge_guard_cells_into_bulk( Field<T, DField, DGrid>& field, const mpi::CartComm& comm );
}

#endif
