#include <mpi.h>
#include <type_traits>
namespace mpi {

  template <typename T_cvref>
  constexpr MPI_Datatype mpi_datatype() noexcept {
    using T = std::remove_cv_t< std::remove_reference_t< T_cvref > >;
    if constexpr ( std::is_same_v<T, char> ) return MPI_CHAR;
    else if ( std::is_same_v<T, short> ) return MPI_SHORT;
    else if ( std::is_same_v<T, int> ) return MPI_INT;
    else if ( std::is_same_v<T, long> ) return MPI_LONG;

    else if ( std::is_same_v<T, unsigned char> ) return MPI_UNSIGNED_CHAR;
    else if ( std::is_same_v<T, unsigned short> ) return MPI_UNSIGNED_SHORT;
    else if ( std::is_same_v<T, unsigned int> ) return MPI_UNSIGNED;
    else if ( std::is_same_v<T, unsigned long> ) return MPI_UNSIGNED_LONG;

    else if ( std::is_same_v<T, float> ) return MPI_FLOAT;
    else if ( std::is_same_v<T, double> ) return MPI_DOUBLE;

    else if ( std::is_same_v<T, long double> ) return MPI_LONG_DOUBLE;
    else if ( std::is_same_v<T, bool> ) return MPI_CXX_BOOL;
    else return MPI_DATATYPE_NULL;
  }

}
