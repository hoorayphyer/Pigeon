#include "parallel/mpi_datatype.hpp"
#include <mpi.h>

namespace mpi {
  template <>
  MPI_Datatype datatype<bool>(bool*) noexcept { return MPI_CXX_BOOL; }
  template <>
  MPI_Datatype datatype<char>(char*) noexcept { return MPI_CHAR; }
  template <>
  MPI_Datatype datatype<short>(short*) noexcept { return MPI_SHORT; }
  template <>
  MPI_Datatype datatype<int>(int*) noexcept { return MPI_INT; }
  template <>
  MPI_Datatype datatype<long>(long*) noexcept { return MPI_LONG; }
  template <>
  MPI_Datatype datatype<long long int>(long long int*) noexcept { return MPI_LONG_LONG_INT; }

  template <>
  MPI_Datatype datatype<unsigned char>(unsigned char*) noexcept { return MPI_UNSIGNED_CHAR; }
  template <>
  MPI_Datatype datatype<unsigned short>(unsigned short*) noexcept { return MPI_UNSIGNED_SHORT; }
  template <>
  MPI_Datatype datatype<unsigned int>(unsigned int*) noexcept { return MPI_UNSIGNED; }
  template <>
  MPI_Datatype datatype<unsigned long>(unsigned long*) noexcept { return MPI_UNSIGNED_LONG; }
  template <>
  MPI_Datatype datatype<unsigned long long int>(unsigned long long int*) noexcept { return MPI_UNSIGNED_LONG_LONG; }

  template <>
  MPI_Datatype datatype<float>(float*) noexcept { return MPI_FLOAT; }
  template <>
  MPI_Datatype datatype<double>(double*) noexcept { return MPI_DOUBLE; }
  template <>
  MPI_Datatype datatype<long double>(long double*) noexcept { return MPI_LONG_DOUBLE; }
}

#include "mpi_cparticle.cpp"
#include <vector>
namespace mpi {
  MPI_Datatype MPI_CPARTICLE = MPI_DATATYPE_NULL;

  MPI_Datatype create_MPI_CPARTICLE () {
    std::vector<cPtc_t> ptcs(2);
    constexpr int numBlocks = 4;
    MPI_Datatype type[numBlocks] = { datatype<real_t>(), datatype<real_t>(),
                                     datatype<typename particle::Specs<real_t>::state_type>(), datatype<char>() };
    int blocklen[numBlocks] = { particle::Specs<real_t>::Dim, particle::Specs<real_t>::Dim, 1, 1 };
    MPI_Aint disp[numBlocks];

    MPI_Get_address( &(ptcs[0].q()), disp );
    MPI_Get_address( &(ptcs[0].p()), disp+1 );
    MPI_Get_address( &(ptcs[0].state()), disp+2 );
    MPI_Get_address( &(ptcs[0].extra()), disp+3 );

    auto base = disp[0];
    for ( int i = 0; i < numBlocks; ++i )
      disp[i] = MPI_Aint_diff( disp[i], base );

    // first create a tmp type
    MPI_Datatype mdt_tmp;
    MPI_Type_create_struct( numBlocks, blocklen, disp, type, &mdt_tmp );
    // adjust in case of mysterious compiler padding
    MPI_Aint sizeofentry;
    MPI_Get_address( ptcs.data() + 1, &sizeofentry );
    sizeofentry = MPI_Aint_diff(sizeofentry, base);

    MPI_Datatype mdt_ptc;
    MPI_Type_create_resized(mdt_tmp, 0, sizeofentry, &mdt_ptc);
    return mdt_ptc;
  }

  template <>
  MPI_Datatype datatype<cPtc_t>(cPtc_t*) noexcept { return MPI_CPARTICLE; }
}
