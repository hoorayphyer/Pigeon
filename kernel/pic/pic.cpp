#include "simulator.hpp"
#include <stdexcept>
#include <time.h>
#include <cstring>
#include "filesys/filesys.hpp"
#include "logger/logger.hpp"
#include "argparser.hpp"

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


// define here to avoid multiple definition
namespace particle{
  map<Properties> properties;
}

int main(int argc, char** argv) {

  auto cli_args = pic::parse_args(argc, argv); // NOTE currently it's not playing any role

  mpi::initialize(argc, argv);

  { // use block to force destruction of potential mpi communicators before mpi::finalize
    pic::this_run_dir = io::init_this_run_dir( pic::datadir_prefix, data_dirname() );
    auto cart_opt = make_cart(pic::dims, pic::periodic);

    // journaling
    fs::mpido( mpi::world, [&](){
                             if ( !fs::exists("journal.txt") ) return;
                             std::ofstream out;
                             out.open("journal.txt", std::ios_base::app);
                             out << "DataDir := " << pic::this_run_dir << std::endl;
                             out.close();
                             // NOTE fs::rename doesn't work on some platforms because of cross-device link.
                             fs::copy_file( "journal.txt", pic::this_run_dir+"/journal.txt" );
                           } );

    fs::mpido( mpi::world, [&](){
                             fs::create_directories(pic::this_run_dir + "/data");
                             fs::create_directories(pic::this_run_dir + "/logs");
                             fs::create_directories(pic::this_run_dir + "/pigeon");
                             fs::copy_file("CMakeLists.txt", pic::this_run_dir + "/pigeon/CMakeLists.txt" );
                             fs::copy_file("pic.hpp", pic::this_run_dir + "/pigeon/pic.hpp" );
                             fs::copy_file("gen.hpp", pic::this_run_dir + "/pigeon/gen.hpp" );
                           } );
    lgr::file.set_filename( pic::this_run_dir + "/logs/rank" + std::to_string(mpi::world.rank()) + ".log" );

    field::set_up<pic::real_t>();
    particle::set_up_properties();
    particle::set_up<pic::real_t>();

    pic::Simulator< pic::DGrid, pic::real_t, particle::Specs, pic::ShapeF, pic::real_j_t, pic::Metric >
      sim( pic::supergrid, cart_opt, pic::guard );

    int init_timestep = sim.load_initial_condition<pic::InitialCondition>();
    sim.set_rng_seed( init_timestep + mpi::world.rank() );

    for ( int ts = init_timestep; ts < init_timestep + pic::total_timesteps; ++ts ) {
      sim.evolve( ts, pic::dt );
    }
    lgr::file.close();
    // barrier everyone to avoid idles calling mpi::finalize() prematurely
    mpi::world.barrier();
  }

  mpi::finalize();
  return 0;
}
