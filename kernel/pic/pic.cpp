#include "pic/simulator.hpp"
#include <stdexcept>
#include <time.h>
#include <cstring>
#include "logger/logger.hpp"
#include "pic/argparser.hpp"
#include "particle/mpi_particle.hpp"

namespace {
  // local directory for storing data symlinks
#ifdef APPARENT_DATA_DIR
  std::string local_data_dir =
    []() {
      std::string str = APPARENT_DATA_DIR;
      fs::remove_slash(str);
      return str;
    }();
#else
  std::string local_data_dir = "Data";
#endif

  // TODOL what if prefix == local_data_dir??
  std::string init_this_run_dir( std::string prefix, std::string dirname ) {
    // use world root time to ensure uniqueness
    std::string this_run_dir;
    if ( mpi::world.rank() == 0 ) {
      prefix = fs::absolute(prefix);
      fs::remove_slash(prefix);
      fs::remove_slash(dirname);

      // in case of running too frequently within a minute, directories with postfixed numbers are created
      if ( fs::exists(prefix + "/" + dirname) ) {
        for ( int n = 1; ; ++n ) {
          if ( !fs::exists(prefix + "/" + dirname + "-" + std::to_string(n)) ) {
            dirname += "-" + std::to_string(n);
            break;
          }
        }
      }
      this_run_dir = prefix + "/" + dirname;

      fs::create_directories(this_run_dir);
      fs::create_directories(local_data_dir);
      fs::create_directory_symlink(this_run_dir, local_data_dir + "/" + dirname);
    }

    if ( mpi::world.size() > 1 ) {
      char buf[200];
      if ( mpi::world.rank() == 0 ) {
        for ( int i = 0; i < this_run_dir.size(); ++i )
          buf[i] = this_run_dir[i];
        buf[this_run_dir.size()] = '\0';
        mpi::world.broadcast(0, buf, 200);
      }
      else {
        mpi::world.broadcast(0, buf, 200);
        this_run_dir = {buf};
      }
    }

    return this_run_dir;

  }
}


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
  return pic::project_name + "-" + subDir;
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

int main(int argc, char** argv) {
  auto cli_args = pic::parse_args(argc, argv);

  std::optional<std::string> resume_dir;
  if ( cli_args.resume_dir ) {
    resume_dir.emplace(fs::absolute(*cli_args.resume_dir));
  }

  pic::ConfFile_t conf;
  if ( cli_args.config_file ) {
    try {
      conf = toml::parse_file(*cli_args.config_file);
    } catch (const toml::parse_error& err) {
      std::cerr
        << "Error parsing file '"sv << *err.source().path
        << "':\n"sv << err.description()
        << "\n  ("sv << err.source().begin << ")"sv
        << std::endl;
      return 1;
    }
  }
  pic::load_configuration(conf);

  if ( cli_args.is_dry_run ) {
    int retcode = 0;
    std::cout << "Dry Run Checks :=" << std::endl;
#if PIC_DEBUG
    std::cout << "\tDebug" << std::endl;
#else
    std::cout << "\tRelease" << std::endl;
#endif
    std::cout << pic::proofread("\t") << std::endl;
    if ( resume_dir  ) {
      if ( !fs::exists(*resume_dir) ) {
        retcode = 1;
        std::cout << "ERROR : Invalid resume directory : " << *resume_dir << std::endl;
      } else if(resume_dir->find("checkpoints/timestep") == std::string::npos ){
        // if the directory exists but is not one of the checkpoints
        retcode = 1;
        std::cout << "ERROR : Invalid resume directory : " << *resume_dir << ". Specify which checkpoint!" << std::endl;
      } else {
        std::cout << "\tResume from : " << *resume_dir << std::endl;
      }
    }

    return retcode;
  }

  mpi::initialize(argc, argv);

  mpi::commit( mpi::Datatype<particle::Particle<pic::real_t, particle::Specs>>{} );

  { // use block to force destruction of potential mpi communicators before mpi::finalize
    pic::this_run_dir = init_this_run_dir( pic::datadir_prefix, data_dirname() );
    auto cart_opt = make_cart(pic::dims, pic::periodic);

    // journaling
    fs::mpido( mpi::world, [&]() {
                             const std::string official_jnl (pic::this_run_dir+"/journal.txt");
                             std::string jnl;
                             if ( cli_args.journal_file ) {
                               // a journal file is specified
                               jnl = fs::absolute(*cli_args.journal_file);
                               if ( !fs::exists(jnl) ) {
                                 std::cout << "Specified journal doesn't exist. Using default journal instead." << std::endl;
                                 jnl = official_jnl;
                               }
                             } else {
                               // if a journal file is not specified, create one
                               jnl = official_jnl;
                             }
                             std::ofstream out;
                             out.open(jnl, std::ios_base::app); // NOTE app creates new file when jnl doesn't exist
#if PIC_DEBUG
                             out << "BuildType := Debug" << std::endl;
#else
                             out << "BuildType := Release" << std::endl;
#endif
                             out << "DataDir := " << pic::this_run_dir << std::endl;
                             if ( resume_dir )
                               out << "Resume := " << *resume_dir << std::endl;
                             out.close();
                             if ( !fs::equivalent(jnl, official_jnl) ) {
                               // NOTE fs::rename doesn't work on some platforms because of cross-device link.
                               fs::copy_file( jnl, official_jnl );
                             }
                           } );

    fs::mpido( mpi::world, [&](){
                             fs::create_directories(pic::this_run_dir + "/data");
                             fs::create_directories(pic::this_run_dir + "/logs");
                             fs::create_directories(pic::this_run_dir + "/pigeon");
                             fs::copy_file("CMakeLists.txt", pic::this_run_dir + "/pigeon/CMakeLists.txt" );
                             fs::copy_file("pic.hpp", pic::this_run_dir + "/pigeon/pic.hpp" );
                             fs::copy_file("pic_impl.hpp", pic::this_run_dir + "/pigeon/pic_impl.hpp" );
                             if ( cli_args.config_file ) {
                               fs::copy_file(*cli_args.config_file, pic::this_run_dir + "/pigeon/conf.toml" );
                             }
                           } );
    lgr::file.set_filename( pic::this_run_dir + "/logs/rank" + std::to_string(mpi::world.rank()) + ".log" );

    auto properties = pic::set_up_particle_properties();

    if (mpi::world.rank() == 0)
      std::cout << "Initializing simulator..." << std::endl;

    pic::Simulator< pic::DGrid, pic::real_t, particle::Specs, pic::ShapeF, pic::real_j_t >
      sim( pic::supergrid, cart_opt, properties );

    if (mpi::world.rank() == 0)
      std::cout << "Loading initial condition..." << std::endl;
    int init_timestep = sim.load_initial_condition( resume_dir );
    sim.set_rng_seed( init_timestep + mpi::world.rank() );

    if (mpi::world.rank() == 0)
      std::cout << "Launch" << std::endl;
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
