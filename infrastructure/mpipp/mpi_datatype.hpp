#ifndef  _MPI_DATATYPE_HPP_
#define  _MPI_DATATYPE_HPP_

#include <mpi.h>

namespace mpi {
  template < typename T >
  struct Datatype {
    static MPI_Datatype type;
    // NOTE for custom datatypes, specialization of Datatype should have the following function
    static void define(MPI_Datatype&);
  };
}

namespace mpi {
  template < typename T >
  constexpr MPI_Datatype datatype(T* = nullptr) noexcept {
    return Datatype<T>::type;
  }

  template <typename T>
  constexpr MPI_Datatype datatype(const T&) noexcept {
    return datatype((T*)0);
  }

  template <typename T>
  constexpr MPI_Datatype datatype(const T*) noexcept {
    return datatype((T*)0);
  }
}

#endif
