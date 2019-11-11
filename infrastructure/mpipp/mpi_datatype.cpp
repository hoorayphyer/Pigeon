#include "mpipp/mpi_datatype.hpp"
#include <mpi.h>

namespace mpi {
  template <>
  MPI_Datatype Datatype<bool>::type = MPI_CXX_BOOL;
  template <>
  MPI_Datatype Datatype<char>::type = MPI_CHAR;
  template <>
  MPI_Datatype Datatype<short>::type = MPI_SHORT;
  template <>
  MPI_Datatype Datatype<int>::type = MPI_INT;
  template <>
  MPI_Datatype Datatype<long>::type = MPI_LONG;
  template <>
  MPI_Datatype Datatype<long long int>::type = MPI_LONG_LONG_INT;

  template <>
  MPI_Datatype Datatype<unsigned char>::type = MPI_UNSIGNED_CHAR;
  template <>
  MPI_Datatype Datatype<unsigned short>::type = MPI_UNSIGNED_SHORT;
  template <>
  MPI_Datatype Datatype<unsigned int>::type = MPI_UNSIGNED;
  template <>
  MPI_Datatype Datatype<unsigned long>::type = MPI_UNSIGNED_LONG;
  template <>
  MPI_Datatype Datatype<unsigned long long int>::type = MPI_UNSIGNED_LONG_LONG;

  template <>
  MPI_Datatype Datatype<float>::type = MPI_FLOAT;
  template <>
  MPI_Datatype Datatype<double>::type = MPI_DOUBLE;
  template <>
  MPI_Datatype Datatype<long double>::type = MPI_LONG_DOUBLE;
}
