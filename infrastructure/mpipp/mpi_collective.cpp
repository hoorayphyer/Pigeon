#include "mpipp/mpi_collective_impl.hpp"
namespace mpi {
  INSTANTIATE_MPI_COLLECTIVE(char);
  INSTANTIATE_MPI_COLLECTIVE(int);
  INSTANTIATE_MPI_COLLECTIVE(long long);
  INSTANTIATE_MPI_COLLECTIVE(unsigned int);
  INSTANTIATE_MPI_COLLECTIVE(std::size_t);
  INSTANTIATE_MPI_COLLECTIVE(float);
  INSTANTIATE_MPI_COLLECTIVE(double);
  INSTANTIATE_MPI_COLLECTIVE(long double);
  template void Collective_Comm<Comm>::barrier(const char*) const;
  template void Collective_Comm<InterComm>::barrier(const char*) const;
}
