#include "io/silo++.cpp"
#include "io/silo_optlist.cpp"
#include "parallel/mpi++.hpp"
#include "catch2/catch.hpp"
#include <unistd.h>

using namespace mpi;
using namespace silo;

SCENARIO("OptList", "[io][silo]") {
  OptList optlist;
  optlist[DBOPT_TIME] = 147.0;
  optlist[DBOPT_DTIME] = 369.0;
  optlist[DBOPT_CYCLE] = 555;
  // DBOptlist* list = optlist;

}

SCENARIO("", "[io][silo]") {
  if ( world.size() > 1 ) {
    std::string filename = "test.silo";
    std::string varname = "silo_data";
    if ( world.rank() == 0 ) {
      auto dbfile = open<Mode::Write>(filename.c_str());
      int silo_data = 147147;
      int dims[1] = {1};
      DBWrite( dbfile, varname.c_str(), &silo_data, dims, 1, DB_INT );
      sleep(1); // have it wait long enough to create the silo file
      close(dbfile);

      int msg = 137;
      world.send(1, 147, &msg, 1 );
    } else if ( world.rank() == 1 ) {
      int msg;
      world.recv(0, 147, &msg, 1 );
      REQUIRE( msg == 137 );
      DBfile* dbfile = DBOpen( filename.c_str(), DB_HDF5, DB_READ );
      int silo_data = 0;
      DBReadVar(dbfile, varname.c_str(), &silo_data);
      DBClose(dbfile);
      REQUIRE( silo_data == 147147 );
    }
  }
}
