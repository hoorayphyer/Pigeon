#ifndef  _FIELD_COMMUNICATION_HPP_
#define  _FIELD_COMMUNICATION_HPP_

namespace mpi { struct Comm; }

namespace field {
  template < typename, int , int > struct Field;

  template < typename T, int DField, int DGrid >
  void sync_guard_cells( Field<T, DField, DGrid>& field, const mpi::Comm& comm );
  // void send_add_guard_cells( auto& field );
}

#endif
