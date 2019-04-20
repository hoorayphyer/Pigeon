#define CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"
#include "parallel/mpi++.hpp"

int main( int argc, char* argv[] ) {
  mpi::initialize();

  int result = Catch::Session().run( argc, argv );

  mpi::finalize();
  return result;
}
