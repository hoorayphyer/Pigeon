#include "all_in_one.hpp"
#include "particle/load_type.hpp"
#include "particle/array.cpp"
#include "parallel/mpi++.hpp"
#include "dye/dynamic_balance.cpp"
#include <unordered_set>

using namespace particle;
SCENARIO("Test calc_new_nprocs", "[dye]") {
  if ( mpi::world.rank() == 0 ) {
    aio::unif_int<int> nens_gen( 1, 100 );
    aio::unif_int<int> ens_size_gen( 1, 100 );
    // aio::gauss_real<double> load_gen( 1000.0, 500.0 );

    int N = 1000;
    while(N--) {
      int nens = nens_gen();
      load_t load_per_ens = 1; // just use a normalized load different than 0
      load_t target_load = 0;

      std::vector<load_t> loads_and_nprocs(2*nens);
      int max_nprocs = 0;
      for ( int i = 0; i < nens; ++i ) {
        loads_and_nprocs[2*i] = load_per_ens;
        auto ens_size = ens_size_gen();
        max_nprocs += ens_size;
        loads_and_nprocs[2*i+1] = ens_size;
      }

      DYNAMIC_SECTION("When all ensembles have same load, and target load is zero, ensembles should have same size and all processes should be recruited " << N ) {
        auto nprocs_new = dye::impl::calc_new_nprocs(loads_and_nprocs, target_load, max_nprocs);
        int total_procs_used = 0;
        for ( int i = 0; i < nprocs_new.size(); ++i ) {
          total_procs_used += nprocs_new[i];
          REQUIRE( ((nprocs_new[i] == max_nprocs/nens) || (nprocs_new[i] == max_nprocs/nens + 1)) );
        }
        REQUIRE( total_procs_used == max_nprocs );
      }

      DYNAMIC_SECTION("When some ensembles have 0 load while others have same load, then there should be at least one proc in each ensemble " << N ) {
        // we skip i == 0 so as to make sure there will be nonzero total load
        for ( int i = 1; i < nens; ++i ) {
          if ( nens_gen() % 4 != 0 ) continue;
          loads_and_nprocs[2*i] = 0;
        }

        std::unordered_set<int> ens_wo_load; // ens_without_load
        for ( int i = 0; i < nens; ++i ) {
          if ( 0 == loads_and_nprocs[2*i] )
            ens_wo_load.insert(i);
        }

        auto nprocs_new = dye::impl::calc_new_nprocs(loads_and_nprocs, target_load, max_nprocs);
        int total_procs_used = 0;
        for ( int i = 0; i < nprocs_new.size(); ++i ) {
          total_procs_used += nprocs_new[i];
          if ( ens_wo_load.find(i) != ens_wo_load.end() ) {
            REQUIRE( nprocs_new[i] == 1 );
          } else {
            auto nprocs_ave = (max_nprocs - ens_wo_load.size());
            if ( ens_wo_load.size() == nens ) nprocs_ave = 0;
            else nprocs_ave /= (nens - ens_wo_load.size());
            CAPTURE(nprocs_ave);
            REQUIRE( (nprocs_new[i] == nprocs_ave || nprocs_new[i] == nprocs_ave + 1) );
          }
        }
        REQUIRE( total_procs_used == max_nprocs );
      }
    }
  }
}

TEMPLATE_TEST_CASE( "Test dynamic balancing on some cartesian topology initially with trivial ensembles","[dye][mpi][.]"
                    , (aio::IndexType<2,1>)
                    ) {
  constexpr int DGrid = 2;
  int nens = 1;
  std::vector<int> cart_dims;
  for ( auto i : TestType::get() ) {
    cart_dims.push_back( i > 0 ? i : -i );
    nens *= cart_dims.back();
  }
  if ( mpi::world.size() > nens ) {
    unsigned int target_load = 0;
    std::vector<load_t> loads(nens);
    // loads = { 1000, 1000, 1000, 1000 };
    loads = { 1000, 1000 };
    auto cart_opt = aio::make_cart( cart_dims, {false,false}, mpi::world );
    auto ens_opt = dye::create_ensemble<DGrid>(cart_opt);
    map<array<double, Specs>> ptcs;
    if ( ens_opt ) {
      auto myload = loads[ens_opt->label()];
      ptcs[species::electron].resize( myload );
    }

    dye::dynamic_load_balance(ptcs, ens_opt, cart_opt, target_load );

  }

}
