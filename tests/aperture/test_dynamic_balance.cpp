#include <iostream>
#include "aperture/dynamic_balance.cpp"
#include "parallel/mpi_datatype.cpp"
#include "parallel/mpi++.cpp"
#include "parallel/mpi_communication.cpp"
#include "catch2/catch.hpp"
#include <random>
#include <ctime>
#include <cmath>
#include <algorithm> // std::min, std::max

using namespace aperture;

SCENARIO("Test get_proc_surplus", "[aperture]") {
  // Most important is to have a feasible surplus
  // std::default_random_engine eng(std::time(nullptr));
  std::default_random_engine eng;
  std::uniform_int_distribution<unsigned int> dist_nproc( 1, 100 );
  std::uniform_real_distribution<double> dist_lg_load( 6.0, 7.0 );

  unsigned int nens = 10;
  load_t target_load = 0;
  unsigned int max_num_procs = 0;

  int N = 10000;

  while( N-- ) {

    // prepare loads
    std::vector<load_t> loads( 2 * nens );
    for ( int i = 0; i < nens; ++i ) {
      loads[2*i] = static_cast<load_t>( std::pow( 10.0, dist_lg_load(eng) ) );
      loads[2*i + 1] = dist_nproc(eng);
    }
    // std::vector<load_t> loads = { 183263, 76, 1162878, 22, 2278921, 68,
    //                               584795, 52, 117258, 6, 2199372, 1,
    //                               136045 , 42, 1506445, 94, 1132029, 10,
    //                               679201 , 71 };
    unsigned int tot_nprocs_old = 0;
    load_t tot_load = 0;
    load_t avld_min_old = 1000*1000*1000*1000; // should be enough
    load_t avld_max_old = 0;
    for ( int i = 0; i < nens; ++i ) {
      tot_load += loads[2*i];
      tot_nprocs_old += loads[2*i + 1];
      auto avld = loads[2*i] / loads[2*i + 1];
      avld_min_old = std::min<load_t>( avld_min_old, avld );
      avld_max_old = std::max<load_t>( avld_max_old, avld );
    }

    { // test
      auto nprocs_new = impl::calc_new_nprocs( loads, target_load, max_num_procs );
      REQUIRE(nprocs_new.size() == nens );
      THEN ("Necessary condition for validity: no nprocs overflow, and each ensemble has at least 1 member") {
        unsigned int tot_nprocs_new = 0;
        for ( auto x : nprocs_new ) {
          REQUIRE(x >= 1);
          tot_nprocs_new += x;
        }
        REQUIRE( tot_nprocs_new <= std::max<unsigned int>(tot_nprocs_old, max_num_procs) );
      }
      load_t avld_min_new = 1000*1000*1000*1000; // should be enough
      load_t avld_max_new = 0;
      for ( int i = 0; i < nens; ++i ) {
        auto avld = loads[2*i] / nprocs_new[i];
        avld_min_new = std::min<load_t>( avld_min_new, avld );
        avld_max_new = std::max<load_t>( avld_max_new, avld );
      }
      CAPTURE(avld_min_new);
      CAPTURE(avld_max_new);
      double load_bal = ((double)avld_max_new - avld_min_new) / ( (double)avld_max_new + avld_min_new );
      REQUIRE( load_bal < 0.3 ); // TODO study how this changes with loads and number of ensembles
    }
  }

}
