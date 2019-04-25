#include "all_in_one.hpp"
#include "io/silo++.hpp"
#include <silo.h>
#include "parallel/mpi++.hpp"
#include "utility/filesys.hpp"
#include <unistd.h>

using namespace mpi;
using namespace silo;
using namespace util;

SCENARIO("create files", "[io][silo][.]") {
  if ( world.rank() == 0 ) {
    file_t dbfile;
    REQUIRE_FALSE(fs::exists("pubg.silo"));
    dbfile = open<Mode::Write>( "pubg.silo" );
    REQUIRE(fs::exists("pubg.silo"));
    fs::remove_all("pubg.silo");
  }
  world.barrier();
}

SCENARIO("OptList", "[io][silo]") {
  if ( world.rank() == 0 ) {
    OptList optlist;
    optlist[DBOPT_TIME] = 147.0f; // float option
    optlist[DBOPT_DTIME] = 369.0; // double option
    optlist[DBOPT_CYCLE] = 555; // int option
    optlist[DBOPT_LO_OFFSET] = std::vector<int>{124,18,384}; // int array
    DBoptlist* list = optlist;

    float* float_opt = (float*)DBGetOption(list, DBOPT_TIME);
    REQUIRE( float_opt[0] == 147.0);

    double* double_opt = (double*)DBGetOption(list, DBOPT_DTIME);
    REQUIRE( double_opt[0] == 369.0);

    int* int_opt = (int*)DBGetOption(list, DBOPT_CYCLE);
    REQUIRE( int_opt[0] == 555);

    int* int_arr_opt = (int*)DBGetOption(list, DBOPT_LO_OFFSET);
    REQUIRE( int_arr_opt[0] == 124 );
    REQUIRE( int_arr_opt[1] == 18 );
    REQUIRE( int_arr_opt[2] == 384 );

  }
  world.barrier();

}

// SCENARIO("", "[io][silo]") {
//   if ( world.size() > 1 && world.rank() < 2 ) {
//     std::string filename = "test.silo";
//     std::string varname = "silo_data";
//     if ( world.rank() == 0 ) {
//       auto dbfile = open<Mode::Write>(filename.c_str());
//       int silo_data = 147147;
//       int dims[1] = {1};
//       DBWrite( dbfile, varname.c_str(), &silo_data, dims, 1, DB_INT );
//       sleep(1); // have it wait long enough to create the silo file
//       close(dbfile);

//       int msg = 137;
//       world.send(1, 147, &msg, 1 );
//     } else if ( world.rank() == 1 ) {
//       int msg;
//       world.recv(0, 147, &msg, 1 );
//       REQUIRE( msg == 137 );
//       DBfile* dbfile = DBOpen( filename.c_str(), DB_HDF5, DB_READ );
//       int silo_data = 0;
//       DBReadVar(dbfile, varname.c_str(), &silo_data);
//       DBClose(dbfile);
//       REQUIRE( silo_data == 147147 );
//     }
//   }

// }
