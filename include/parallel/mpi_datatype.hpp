#ifndef  _MPI_DATATYPE_HPP_
#define  _MPI_DATATYPE_HPP_

typedef struct ompi_datatype_t *MPI_Datatype;

namespace mpi {
  template <typename Type>
  MPI_Datatype datatype(Type* = nullptr) noexcept;

  template <typename T>
  MPI_Datatype datatype(const T&) noexcept {
    return datatype((T*)0);
  }
}

#endif
