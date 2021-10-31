#pragma once

#include <vector>

#include "mpipp/mpi_datatype.hpp"
#include "particle/particle.hpp"

namespace mpi {
template <typename T, template <typename> class S>
struct Datatype<particle::Particle<T, S>> {
  static MPI_Datatype type;

  static void define(MPI_Datatype& dtype) {
    using namespace particle;
    std::vector<Particle<T, S>> ptcs(2);
    constexpr int numBlocks = 4;
    MPI_Datatype type[numBlocks] = {Datatype<T>::type, Datatype<T>::type,
                                    Datatype<T>::type,
                                    Datatype<typename S<T>::state_type>::type};
    int blocklen[numBlocks] = {S<T>::Dim, S<T>::Dim, 1, 1};
    MPI_Aint disp[numBlocks];

    MPI_Get_address(&(ptcs[0].q()), disp);
    MPI_Get_address(&(ptcs[0].p()), disp + 1);
    MPI_Get_address(&(ptcs[0].frac()), disp + 2);
    MPI_Get_address(&(ptcs[0].state()), disp + 3);

    auto base = disp[0];
    for (int i = 0; i < numBlocks; ++i) disp[i] = MPI_Aint_diff(disp[i], base);

    // first create a tmp type
    MPI_Datatype mdt_tmp;
    MPI_Type_create_struct(numBlocks, blocklen, disp, type, &mdt_tmp);
    // adjust in case of mysterious compiler padding
    MPI_Aint sizeofentry;
    MPI_Get_address(ptcs.data() + 1, &sizeofentry);
    sizeofentry = MPI_Aint_diff(sizeofentry, base);

    MPI_Type_create_resized(mdt_tmp, 0, sizeofentry, &dtype);

    MPI_Type_free(&mdt_tmp);
  }
};

template <typename T, template <typename> class S>
MPI_Datatype Datatype<particle::Particle<T, S>>::type = MPI_DATATYPE_NULL;
}  // namespace mpi
