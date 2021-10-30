#include "mpipp/mpi_p2p_impl.hpp"
namespace mpi {
INSTANTIATE_MPI_P2P(char);
INSTANTIATE_MPI_P2P(int);
INSTANTIATE_MPI_P2P(unsigned int);
INSTANTIATE_MPI_P2P(float);
INSTANTIATE_MPI_P2P(double);
INSTANTIATE_MPI_P2P(long double);
}  // namespace mpi
