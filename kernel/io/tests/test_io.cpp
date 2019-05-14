#include "testfw/testfw.hpp"
#include "io/io_impl.hpp"
#include <filesystem>

using namespace std::filesystem;
using namespace io;

SCENARIO( "Test init data export directory", "[io][mpi]" ) {
  std::string prefix = "TDETDE/REAL_DATA";
  io::local_data_dir = "TDETDE/LOCAL";
  if ( mpi::world.rank() == 0 ) {
    remove_all(prefix);
    remove_all(local_data_dir);
  }
  // NOTE dirname is significant only on rank 0
  std::string dirname = "FOO-M416-" + std::to_string(mpi::world.rank());
  auto this_run_dir = init_this_run_dir( prefix , dirname );

  REQUIRE( equivalent(this_run_dir, prefix + "/FOO-M416-0") );

  std::vector<path> paths;
  paths.emplace_back(prefix);
  REQUIRE( exists(paths.back()) );
  REQUIRE( is_directory(paths.back()) );

  paths.emplace_back( prefix + "/" + dirname);
  if ( mpi::world.rank() == 0 ) {
    REQUIRE( exists(paths.back()) );
    REQUIRE( is_directory(paths.back()) );
  } else {
    REQUIRE_FALSE( exists(paths.back()) );
  }

  paths.emplace_back(io::local_data_dir);
  REQUIRE( exists(paths.back()) );
  REQUIRE( is_directory(paths.back()) );

  paths.emplace_back(io::local_data_dir + "/" + dirname);
  if ( mpi::world.rank() == 0 ) {
    REQUIRE( exists(paths.back()) );
    REQUIRE( is_symlink(paths.back()) );
  } else {
    REQUIRE_FALSE( exists(paths.back()) );
  }
}

SCENARIO( "Test downsample", "[io]") {
  if ( mpi::world.rank() == 0 ) {

  }
}
