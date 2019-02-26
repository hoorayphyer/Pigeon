#include "parallel/mpi++.cpp"
#include "parallel/mpi_communication.cpp"
#include "catch2/catch.hpp"

using namespace mpi;

SCENARIO("World", "[parallel][mpi]") {
  if ( world.size() == 2 ) {
    int myrank = world.rank();
    // std::cout << "world size = " << world.size() << " ," << myrank << std::endl;
    if ( 0 == myrank ) {
      int msg = 147;
      world.send( 1, 22, &msg );
    } else {
      int msg;
      world.recv( 0, 22, &msg, 1 );
      REQUIRE( msg == 147 );
    }
  }
}

SCENARIO("Group", "[parallel][mpi]") {
  auto grp_world = world.group();
  SECTION("Constructors") {
    WHEN("Default") {
      Group grp;
      REQUIRE( MPI_UNDEFINED == grp.rank() );
      REQUIRE( 0 == grp.size() );
    }

    WHEN("From MPI_WORLD_GROUP") {
      AND_WHEN("no duplicates") {
      }
      AND_WHEN("there are duplicates") {
      }
    }
  }

  SECTION("Operations") {

  }

}

SCENARIO("Comm", "[parallel][mpi]") {
  Comm comm( world, world.group() );
  if ( world.size() == 2 ) {
    int myrank = comm.rank();
    if ( 0 == myrank ) {
      int msg = 147;
      comm.send( 1, 22, &msg );
    } else {
      int msg;
      comm.recv( 0, 22, &msg, 1 );
      REQUIRE( 147 == msg );
    }
  }
}

// SCENARIO("intra P2P communications", "[parallel][mpi]") {}

// SCENARIO("intra collective communications", "[parallel][mpi]") {}

SCENARIO("intercomm P2P communications", "[parallel][mpi]") {
  if ( world.size() == 2 ) {
    SECTION("each comm has one member") {
      int myrank = world.rank();
      Comm local ( world, myrank, myrank );
      std::optional<Comm> peer_comm(world);
      InterComm inter( local, 0, peer_comm, 1-myrank, 147 );
      if ( 0 == myrank ) {
        int msg = 147;
        inter.send( 0, 22, &msg, 1 );
      } else {
        int msg;
        inter.recv( 0, 22, &msg, 1 );
        REQUIRE(msg == 147);
      }
    }
  }

}
