#include "simulator.hpp"
#include <stdexcept>
#include <time.h>
#include <cstring>
#include "filesys/filesys.hpp"
#include "logger/logger.hpp"

std::string data_dirname() {
  char subDir[100] = {};
  for ( int i = 0; i < 100; ++i )
    subDir[i] = '\0';
  if ( mpi::world.rank() == 0 ) {
    char myTime[100] = {};
    time_t rawtime;
    struct tm* timeinfo;
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(myTime, 100, "%Y%m%d-%H%M", timeinfo);
    snprintf(subDir, sizeof(subDir), "%s", myTime);
  }
  return std::string(pic::project_name) + "-" + subDir;
}

std::optional<mpi::CartComm> make_cart( const apt::array<int,pic::DGrid>& dims, const apt::array<bool,pic::DGrid> periodic ) {
  std::optional<mpi::CartComm> cart_opt;
  int nprmy = 1;
  for ( auto i : dims )
    nprmy *= i;

  if ( nprmy > mpi::world.size() )
    throw std::invalid_argument("Size of the Cartesian topology exceeds the size of world!");

  // simply use first few members in mpi::world as primaries
  bool is_prmy = mpi::world.rank() < nprmy;
  auto comm_tmp = mpi::world.split( {is_prmy}, mpi::world.rank() );
  if ( is_prmy ) {
    std::vector<int> d(pic::DGrid);
    std::vector<bool> p(pic::DGrid);
    for ( int i = 0; i < pic::DGrid; ++i ) {
      d[i] = dims[i];
      p[i] = periodic[i];
    }
    cart_opt.emplace( *comm_tmp, d, p );
  }
  return cart_opt;
}

int main() {

  mpi::initialize();

  { // use block to force destruction of potential mpi communicators before mpi::finalize
    pic::this_run_dir = io::init_this_run_dir( pic::datadir_prefix, data_dirname() );
    auto cart_opt = make_cart(pic::dims, pic::periodic);

    fs::mpido( mpi::world, [&](){
                             fs::create_directories(pic::this_run_dir + "/data");
                             fs::create_directories(pic::this_run_dir + "/logs");
                           } );
    lgr::file.open( pic::this_run_dir + "/logs/rank" + std::to_string(mpi::world.rank()) + ".log" );

    pic::Simulator< pic::DGrid, pic::real_t, particle::Specs, pic::ShapeF, pic::real_j_t, pic::Metric >
      sim( pic::supergrid, cart_opt, pic::guard );

    int init_timestep = sim.load_initial_condition<pic::InitialCondition>();
    sim.set_rng_seed( init_timestep + mpi::world.rank() );

    for ( int ts = init_timestep; ts < init_timestep + pic::total_timesteps; ++ts ) {
      lgr::file % "==== Timestep " << ts << " ====" << std::endl;
      lgr::file.indent_append("\t");
      sim.evolve( ts, pic::dt );
      lgr::file.indent_reset();

      // occasionally barrier everyone to avoid idles running too fast
      if ( ts % 100 == 0 ) mpi::world.barrier();
    }
    lgr::file.close();
    mpi::world.barrier();
  }

  mpi::finalize();
  return 0;
}
