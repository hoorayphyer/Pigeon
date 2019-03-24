#include "aperture.hpp"
#include "traits.hpp"
#include <stdexcept>

using namespace aperture;
using namespace traits;

// TODO
// namespace mpi { // Define MPI_CPARTICLE
//   template < typename T, int DPtc, typename state_t >
//   MPI_Datatype datatype_for_cParticle() {
//     std::vector< particle::cParticle< T, DPtc, state_t > ptcs(2);
//     constexpr int numBlocks = 3;
//     MPI_Datatype type[numBlocks] = { datatype<T>(), datatype<T>(), datatype<state_t>() };
//     int blocklen[numBlocks] = { DPtc, DPtc, 1 };
//     MPI_Aint disp[numBlocks];

//     MPI_Get_address( &(ptcs[0].q()), disp );
//     MPI_Get_address( &(ptcs[0].p()), disp+1 );
//     MPI_Get_address( &(ptcs[0].state()), disp+2 );

//     auto base = disp[0];
//     for ( int i = 0; i < numBlocks; ++i )
//       disp[i] = MPI_Aint_diff( disp[i], base );

//     // first create a tmp type
//     MPI_Datatype mdt_tmp;
//     MPI_Type_create_struct( numBlocks, blocklen, disp, type, &mdt_tmp );
//     // adjust in case of mysterious compiler padding
//     MPI_Aint sizeofentry;
//     MPI_Get_address( ptcs.data() + 1, &sizeofentry );
//     sizeofentry = MPI_Aint_diff(sizeofentry, base);

//     MPI_Datatype mdt_ptc;
//     MPI_Type_create_resized(mdt_tmp, 0, sizeofentry, &mdt_ptc);
//     return mdt_ptc;
//   }

//   MPI_Datatype MPI_CPARTICLE = datatype_for_cParticle<cParticle<T,DPtc,state_t>>();

//   // MPI_Type_commit(&mpi::MPI_CPARTICLE);
//   // MPI_Type_free( &mpi::MPI_CPARTICLE );
// }

constexpr auto make_supergrid() {
  knl::Grid< real_t, DGrid > grid;
  apt::foreach<0,DGrid>
    ( []( auto& g, auto l, auto d ) {
        g = {l[LFT], l[RGT], d};
      }, grid, grid_limits, grid_dims );
  return grid;
}

int main() {

  mpi::initialize();
  // TODO ptc updater constructor

  std::vector<int> dims (DGrid);
  std::vector<bool> periodic (DGrid);
  int total_timesteps = 10;
  real_t dt = 0.001;
  real_t unit_e = 100;

  std::optional<mpi::CartComm> cart_opt;
  // TODO TODO
  // {
  //   int nprmy = 1;
  //   for ( auto i : dims )
  //     nprmy *= i;

  //   if ( nprmy > mpi::world.size() )
  //     throw std::invalid_argument("Size of the Cartesian topology exceeds the size of world!");

  //   // simply use first few members in mpi::world as primaries
  //   bool is_prmy = mpi::world.rank() < nprmy;
  //   auto comm_tmp = mpi::world.split( {is_prmy} );
  //   if ( is_prmy ) cart_opt.emplace( *comm_tmp, dims, periodic );
  // }

  constexpr auto supergrid = make_supergrid();

  using Aperture = Aperture< real_t, DGrid, ptc_state_t >;
  Aperture aperture( cart_opt, supergrid, guard );
  // TODO initial condition. Hard Code for now. Will set the following
  int timestep_begin = 0;

  for ( int ts = timestep_begin; ts < timestep_begin + total_timesteps; ++ts ) {
    aperture.evolve( ts, dt, unit_e );
  }

  mpi::finalize();
  return 0;
}
