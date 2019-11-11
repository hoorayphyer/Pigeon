#define CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"
#include "mpipp/mpi++.hpp"

int main( int argc, char* argv[] ) {
  mpi::initialize(argc, argv);

  int result = 0;
  {
    result = Catch::Session().run( argc, argv );
  }

  mpi::finalize();
  return result;
}
