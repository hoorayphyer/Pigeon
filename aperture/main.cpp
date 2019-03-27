#include "aperture.hpp"
#include "traits.hpp"
#include <stdexcept>

using namespace aperture;
using namespace traits;

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

  std::vector<int> dims (DGrid);
  std::vector<bool> periodic (DGrid);
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

  constexpr auto supergrid = make_supergrid();

  // TODO initial condition. Hard Code for now. Will set the following
  int timestep_begin = 0;

  util::Rng<real_t> rng{};
  rng.set_seed( timestep_begin + mpi::world.rank() );

  using Aperture = Aperture< real_t, DGrid, ptc_state_t >;
  Aperture aperture( supergrid, cart_opt, guard, rng );

  for ( int ts = timestep_begin; ts < timestep_begin + total_timesteps; ++ts ) {
    aperture.evolve( ts, dt, unit_e );
  }

  mpi::finalize();
  return 0;
}
