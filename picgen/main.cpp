#include "picgen.hpp"
#include <stdexcept>

using namespace pic;

int main() {

  mpi::initialize();

  { // use block to force destruction of potential mpi communicators before mpi::finalize
    std::vector<int> dims = { 1, 1 };
    std::vector<bool> periodic = {false,false};
    int total_timesteps = 10;
    real_t dt = 0.001;
    real_t unit_e = 100;

    std::optional<mpi::CartComm> cart_opt;
    {
      int nprmy = 1;
      for ( auto i : dims )
        nprmy *= i;

      if ( nprmy > mpi::world.size() )
        throw std::invalid_argument("Size of the Cartesian topology exceeds the size of world!");

      // simply use first few members in mpi::world as primaries
      bool is_prmy = mpi::world.rank() < nprmy;
      auto comm_tmp = mpi::world.split( {is_prmy}, mpi::world.rank() );
      if ( is_prmy ) cart_opt.emplace( *comm_tmp, dims, periodic );
    }

    constexpr knl::Grid<real_t,DGrid> supergrid
      = {{ { 0.0, std::log(30.0), 16 }, { 0.0, PI, 16 } }};
    constexpr int guard = 1;

    // TODO initial condition. Hard Code for now. Will set the following
    int timestep_begin = 0;

    util::Rng<real_t> rng{};
    rng.set_seed( timestep_begin + mpi::world.rank() );

    using PIC = PIC< DGrid, real_t, particle::Specs, ShapeF, real_j_t, coordinate_system >;
    PIC simulation( supergrid, cart_opt, guard, rng, unit_e );

    for ( int ts = timestep_begin; ts < timestep_begin + total_timesteps; ++ts ) {
      simulation.evolve( ts, dt );
    }
  }

  mpi::finalize();
  return 0;
}
