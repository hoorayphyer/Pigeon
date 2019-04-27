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
    REQUIRE_FALSE(fs::exists("pubg__cc.silo"));
    dbfile = open<Mode::Write>( "pubg__cc.silo" );
    REQUIRE(fs::exists("pubg__cc.silo"));
    fs::remove_all("pubg__cc.silo");
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

SCENARIO("putters", "[io][silo][.]") {
  if ( world.rank() == 0 ) {
    WHEN("put a quadmesh into a silo file") {
      auto dbfile = open<Mode::Write>( "pubg.silo" );

      apt::array<int,3> dims = { 16, 32, 64 };
      std::vector<std::vector<double>> coords(dims.size());
      for ( int i = 0; i < dims.size(); ++i ) {
        coords[i].resize(dims[i]);
        for ( auto& x : coords[i] ) x = dims[i] / 3.0;
      }

      OptList optlist;
      optlist[DBOPT_TIME] = 147.0f; // float option
      optlist[DBOPT_DTIME] = 369.0; // double option
      optlist[DBOPT_CYCLE] = 555; // int option
      optlist[DBOPT_LO_OFFSET] = std::vector<int>{4,10,2}; // int array
      optlist[DBOPT_HI_OFFSET] = std::vector<int>{1,0,2}; // int array

      dbfile.put_mesh( "mesh", coords, optlist );
      silo::close(dbfile);

      AND_WHEN("later reopened for read") {
        auto dbfile = open<Mode::Read>( "pubg.silo" );
        auto* quadmesh = DBGetQuadmesh(dbfile, "mesh");
        REQUIRE( quadmesh->ndims == 3 );
        REQUIRE( quadmesh->datatype == DB_DOUBLE );
        REQUIRE( quadmesh->dims[0] == 16 );
        REQUIRE( quadmesh->dims[1] == 32 );
        REQUIRE( quadmesh->dims[2] == 64 );
        REQUIRE( quadmesh->coordtype == DB_COLLINEAR );
        double* coord;
        coord = (double*)quadmesh->coords[0];
        for ( int i = 0; i < 16; ++i ) REQUIRE( coord[0] == 16 / 3.0 );
        coord = (double*)quadmesh->coords[1];
        for ( int i = 0; i < 32; ++i ) REQUIRE( coord[1] == 32 / 3.0 );
        coord = (double*)quadmesh->coords[2];
        for ( int i = 0; i < 64; ++i ) REQUIRE( coord[2] == 64 / 3.0 );

        // check optlist
        REQUIRE( quadmesh->cycle == 555 );
        REQUIRE( quadmesh->time == 147.0f );
        REQUIRE( quadmesh->dtime == 369.0 );
        // these come from LO and HI offsets
        REQUIRE( quadmesh->min_index[0] == 4 );
        REQUIRE( quadmesh->max_index[0] == 16 - 1 - 1 );
        REQUIRE( quadmesh->min_index[1] == 10 );
        REQUIRE( quadmesh->max_index[1] == 32 - 1 - 0 );
        REQUIRE( quadmesh->min_index[2] == 2 );
        REQUIRE( quadmesh->max_index[2] == 64 - 1 - 2 );

        silo::close(dbfile);
      }

      fs::remove_all("pubg.silo");
    }
  }
  world.barrier();
}

SCENARIO("pmpio putters", "[io][silo]") {
  if ( world.rank() == 0 ) {
    WHEN("put a quadmesh into a silo file") {
      auto dbfile = pmpio::open<Mode::Write>( "pmpio_pubg.silo", "test_pmpio", world, 1 );

      apt::array<int,3> dims = { 16, 32, 64 };
      std::vector<std::vector<double>> coords(dims.size());
      for ( int i = 0; i < dims.size(); ++i ) {
        coords[i].resize(dims[i]);
        for ( auto& x : coords[i] ) x = dims[i] / 3.0;
      }

      OptList optlist;
      optlist[DBOPT_TIME] = 147.0f; // float option
      optlist[DBOPT_DTIME] = 369.0; // double option
      optlist[DBOPT_CYCLE] = 555; // int option
      optlist[DBOPT_LO_OFFSET] = std::vector<int>{4,10,2}; // int array
      optlist[DBOPT_HI_OFFSET] = std::vector<int>{1,0,2}; // int array

      dbfile.put_mesh( "mesh", coords, optlist );
      silo::close(dbfile);

      AND_WHEN("later reopened for read") {
        auto dbfile = open<Mode::Read>( "pmpio_pubg.silo" );
        auto* quadmesh = DBGetQuadmesh(dbfile, "test_pmpio/mesh");
        REQUIRE( quadmesh->ndims == 3 );
        REQUIRE( quadmesh->datatype == DB_DOUBLE );
        REQUIRE( quadmesh->dims[0] == 16 );
        REQUIRE( quadmesh->dims[1] == 32 );
        REQUIRE( quadmesh->dims[2] == 64 );
        REQUIRE( quadmesh->coordtype == DB_COLLINEAR );
        double* coord;
        coord = (double*)quadmesh->coords[0];
        for ( int i = 0; i < 16; ++i ) REQUIRE( coord[0] == 16 / 3.0 );
        coord = (double*)quadmesh->coords[1];
        for ( int i = 0; i < 32; ++i ) REQUIRE( coord[1] == 32 / 3.0 );
        coord = (double*)quadmesh->coords[2];
        for ( int i = 0; i < 64; ++i ) REQUIRE( coord[2] == 64 / 3.0 );

        // check optlist
        REQUIRE( quadmesh->cycle == 555 );
        REQUIRE( quadmesh->time == 147.0f );
        REQUIRE( quadmesh->dtime == 369.0 );
        // these come from LO and HI offsets
        REQUIRE( quadmesh->min_index[0] == 4 );
        REQUIRE( quadmesh->max_index[0] == 16 - 1 - 1 );
        REQUIRE( quadmesh->min_index[1] == 10 );
        REQUIRE( quadmesh->max_index[1] == 32 - 1 - 0 );
        REQUIRE( quadmesh->min_index[2] == 2 );
        REQUIRE( quadmesh->max_index[2] == 64 - 1 - 2 );

        silo::close(dbfile);
      }

      fs::remove_all("pmpio_pubg.silo");
    }
  }
  world.barrier();
}

SCENARIO("multiputters with namescheme", "[io][silo]") {
  if ( world.rank() == 0 ) {
    WHEN("put a 3 blocks by 4 blocks into 2 silo files each of which contains 6 directories") {
      const std::vector<int> dims = { 3, 4 };

      const std::string prefix = "pubg_namescheme/";
      fs::create_directories(prefix);

      // create inidividual quadmesh
      for ( int y = 0; y < 4; ++y ) {
        for ( int x = 0; x < 3; ++x ) {
          auto db = open<Mode::Write>( prefix + "set" + std::to_string((x+3*y)%2) );
          std::string dname = "cart_00" + std::to_string(x) + "_00" + std::to_string(y);
          // create a 1x1 mesh
          std::vector<std::vector<double>> coords = {{x},{y}};
          db.put_mesh("submesh", coords );
        }
      }

      std::vector<int> strides ( dims.size() + 1 );
      strides[0] = 1;
      for ( int i = 0; i < dims.size(); ++i ) strides[i+1] = strides[i] * dims[i];

      // set up file_ns and block_ns
      constexpr char delimiter = '|';
      std::string file_ns = delimiter + prefix + "set%d.silo" + delimiter + "n%2";
      std::string block_ns = delimiter + std::string("cart_%03d_%03d");

      for ( int i = 0; i < dims.size(); ++i ) {
        block_ns += delimiter + std::string("(n%" + std::to_string(strides[i+1]) + ")/") + std::to_string(strides[i]);
      }

      auto dbfile = open<Mode::Write>( prefix + "pubg.silo" );
      dbfile.put_multimesh( "multimesh", strides.back(), file_ns, block_ns );

      fs::remove_all(prefix);
    }
  }
  world.barrier();
}
