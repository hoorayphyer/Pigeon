#include "testfw/testfw.hpp"
#include "particle/load_type.hpp"
#include "particle/array_impl.hpp"
#include "mpipp/mpi++.hpp"
#include "dye/ensemble_impl.hpp"
#include "dye/dynamic_balance_impl.hpp"
#include "pic.hpp"
#include <unordered_set>

using namespace particle;
SCENARIO("Test calc_new_nprocs", "[dye][.]") {
  if ( mpi::world.rank() == 0 ) {
    aio::unif_int<int> nens_gen( 1, 100 );
    aio::unif_int<int> ens_size_gen( 1, 100 );

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

TEST_CASE("Test calc_new_nprocs with nonzero target load", "[dye][.]") {
  if ( mpi::world.rank() == 0 ) {
    std::vector<load_t> loads_and_nprocs = { 584,1,588,1,584,1,588,1,
                                             0,1,0,1, 0,1,0,1,
                                             0,1,0,1, 0,1,0,1,
                                             0,1,0,1,0,1,0,1 };
    load_t target_load = 100000;
    int max_nprocs = 28;
    auto nprocs_new = dye::impl::calc_new_nprocs(loads_and_nprocs, target_load, max_nprocs);
    for ( auto x : nprocs_new )
      REQUIRE(x == 1);
  }
}

TEMPLATE_TEST_CASE( "Test bifurcate","[dye][mpi][.]"
                    , (std::integral_constant<int,4>)
                    , (std::integral_constant<int,7>)
                    ) {
  constexpr auto nprocs = TestType::value;
  if ( nprocs > 1 && mpi::world.size() >= nprocs ) {
    auto parent = *(mpi::world.split( mpi::world.rank() < nprocs ));
    // what bifurcate does is it splits an intra comm into two separate intra comms and creates an intercomm between them
    WHEN("all members don't share the same color") {
      if ( mpi::world.rank() < nprocs ) {
        bool color = parent.rank() < nprocs / 2;
        auto[ intra, inter ] = dye::impl::bifurcate( parent, color );
        REQUIRE(inter);
        if ( color ) {
          REQUIRE( intra.size() == nprocs / 2 );
          REQUIRE( inter->remote_size() == parent.size() - nprocs / 2 );
        } else {
          REQUIRE( intra.size() == parent.size() - nprocs / 2 );
          REQUIRE( inter->remote_size() == nprocs / 2 );
        }
      }
      mpi::world.barrier();
    }

    WHEN("all members share same color") {
      if ( mpi::world.rank() < nprocs ) {
        bool color = true;
        auto[ intra, inter ] = dye::impl::bifurcate( parent, color );
        REQUIRE_FALSE(inter);
        REQUIRE(intra.size() == parent.size());
      }
      mpi::world.barrier();
    }
  }
}

TEMPLATE_TEST_CASE( "Test relinguish_data","[dye][mpi][.]"
                    , (std::integral_constant<int,1>)
                    , (std::integral_constant<int,2>)
                    , (std::integral_constant<int,4>)
                    , (std::integral_constant<int,8>)
                    ) {
  // the sending side of an intercommunicator consists of MPI_ROOT and one MPI_PROC_NULL
  // the receiving side consists of number of processes specified by the template parameter
  constexpr auto nremotes = TestType::value;
  static_assert(nremotes > 0);
  if ( mpi::world.size() >= nremotes + 2 ) {
    auto parent = *( mpi::world.split( mpi::world.rank() < nremotes + 2 ) );
    if ( mpi::world.rank() < nremotes + 2 ) {
      auto[intra, itc] = dye::impl::bifurcate( parent, mpi::world.rank() < 2 );
      array<double,Specs> ptcs;
      if ( mpi::world.rank() == 0 ) {
        ptcs.resize(1000);
        for ( int i = 0; i < ptcs.size(); ++i ) {
          ptcs[i].q()[0] = 132.0;
          ptcs[i].p()[0] = -546.0;
          ptcs[i].set(flag::secondary);
        }
        dye::impl::relinguish_data( ptcs, *itc, MPI_ROOT );
        REQUIRE( ptcs.size() == 0 );
      } else if ( mpi::world.rank() == 1 ) {
        ptcs.resize(10);
        dye::impl::relinguish_data( ptcs, *itc, MPI_PROC_NULL );
        REQUIRE( ptcs.size() == 10 );
      } else {
        dye::impl::relinguish_data( ptcs, *itc, 0 );
        REQUIRE( ptcs.size() == 1000 / intra.size() );
        for ( int i = 0; i < ptcs.size(); ++i ) {
          REQUIRE( ptcs[i].q()[0] == 132.0 );
          REQUIRE( ptcs[i].p()[0] == -546.0 );
          REQUIRE( ptcs[i].is(flag::secondary) );
        }
      }
    }
  }
}

TEMPLATE_TEST_CASE( "Test assign_labels between primaries and idles","[dye][mpi][.]"
                    , (std::integral_constant<int,4>)
                    , (std::integral_constant<int,7>)
                    ) {
  constexpr auto nprocs = TestType::value;
  if ( nprocs > 1 && mpi::world.size() >= nprocs ) {
    auto prmy_idle_comm = *(mpi::world.split( mpi::world.rank() < nprocs ));
    if ( mpi::world.rank() < nprocs ) {
      // what bifurcate does is it splits an intra comm into two separate intra comms and creates an intercomm between them
      int nprmy = nprocs / 2;
      bool is_primary = prmy_idle_comm.rank() < nprmy;
      std::vector<int> deficits; // significant only at primaries
      if ( is_primary ) {
        deficits = std::vector<int>(nprmy,0);
        for ( int i = 0; i < nprocs - nprmy; ++i ) ++deficits[i % nprmy];
      }

      std::optional<int> cur_label; // current label
      if ( is_primary ) cur_label.emplace( prmy_idle_comm.rank() );

      auto[ intra, job_market ] = dye::impl::bifurcate( prmy_idle_comm, is_primary );
      REQUIRE(job_market);

      auto new_label = dye::impl::assign_labels( job_market, deficits, cur_label );
      REQUIRE(new_label);
      if ( is_primary ) REQUIRE( *new_label == *cur_label );
    }
  }
  mpi::world.barrier();
}

TEST_CASE( "Test assign_labels: a specific example ","[dye][mpi][.]") {
  constexpr auto nprocs = 28;
  if ( nprocs > 1 && mpi::world.size() >= nprocs ) {
    auto prmy_idle_comm = *(mpi::world.split( mpi::world.rank() < nprocs ));
    if ( mpi::world.rank() < nprocs ) {
      int nprmy = 16;
      bool is_primary = prmy_idle_comm.rank() < nprmy;
      std::vector<int> deficits; // significant only at primaries
      if ( is_primary ) {
        deficits = std::vector<int>(nprmy,0);
        deficits[0] = 2;
        deficits[1] = 3;
        deficits[2] = 2;
        deficits[3] = 3;
      }

      std::optional<int> cur_label; // current label
      if ( is_primary ) cur_label.emplace( prmy_idle_comm.rank() );

      auto[ intra, job_market ] = dye::impl::bifurcate( prmy_idle_comm, is_primary );
      REQUIRE(job_market);

      auto new_label = dye::impl::assign_labels( job_market, deficits, cur_label );
      if ( is_primary ) REQUIRE( *new_label == *cur_label );

      auto new_intra = mpi::world.split(new_label);
      if (new_intra) new_intra -> barrier();
    }
  }
  mpi::world.barrier();
}

TEMPLATE_TEST_CASE( "Test detailed balance","[dye][mpi][.]"
                    , (std::integral_constant<int,2>)
                    , (std::integral_constant<int,4>)
                    , (std::integral_constant<int,16>)
                    ) {
  constexpr auto ens_size = TestType::value;
  if ( mpi::world.size() >= ens_size ) {
    auto intra = mpi::world.split( mpi::world.rank() < ens_size );
    if ( mpi::world.rank() < ens_size ) {
      aio::gauss_real<double> load_gen( 1000.0, 1000.0 );
      int N = 100;
      while ( N-- ) {
        std::vector<load_t> loads(ens_size);
        if ( intra->rank() == 0 ) {
          for ( auto& l : loads ) {
            auto tmp = load_gen();
            tmp = std::max<double>(tmp, 0);
            tmp = std::min<double>(tmp, 10000);
            l = static_cast<load_t>(tmp);
          }
        }
        intra->broadcast(0,loads.data(), loads.size());

        array<double, Specs> ptcs;
        ptcs.resize(loads[intra->rank()]);

        dye::detailed_balance( ptcs, *intra );

        load_t total_load = 0;
        for ( auto l : loads ) total_load += l;
        CAPTURE( loads, ptcs.size(), total_load / ens_size );
        REQUIRE( (ptcs.size() == total_load / ens_size || ptcs.size() == 1 + total_load / ens_size ) );
      }
    }
    mpi::world.barrier();
  }
}

TEMPLATE_TEST_CASE( "Test dynamic balancing on some cartesian topology initially with trivial ensembles","[dye][mpi][.]"
                    , (aio::IndexType<2,1>)
                    ) {
  // TODO test a shrinking ensemble
  constexpr int DGrid = 2;
  int nens = 1;
  std::vector<int> cart_dims;
  for ( auto i : TestType::get() ) {
    cart_dims.push_back( i > 0 ? i : -i );
    nens *= cart_dims.back();
  }
  // NOTE TODOL current implementation requires explicit touch-create before detailed balance
  map<array<double, Specs>> ptcs;
  map<std::vector<load_t>> loads;
  {
    auto& x = ptcs[species::electron];
    auto& y = ptcs[species::ion];
    auto& xx = loads[species::electron];
    auto& yy = loads[species::ion];
  }


  if ( mpi::world.size() > nens ) {
    aio::gauss_real<double> load_gen( 1000.0, 500.0 );
    unsigned int target_load = 0;
    auto cart_opt = aio::make_cart( cart_dims, {false,false}, mpi::world );
    auto ens_opt = dye::create_ensemble<DGrid>(cart_opt);

    for ( auto&[sp, load] : loads ) {
      load.resize(nens);
      for ( int label = 0; label < nens; ++label ) {
        auto x = load_gen();
        x = std::max<double>(x, 0.0);
        x = std::min<double>(x, 100000.0);
        load[label] = x;

        if ( ens_opt && ens_opt->label() == label )
          ptcs[sp].resize( load[label] );
      }
    }

    dye::dynamic_load_balance(ptcs, ens_opt, cart_opt, target_load );

  }

}

SCENARIO("Stress test", "[dye][mpi][.]") {
  const int num_ens = 1;
  const int num_procs = mpi::world.size();
  const int num_cycles = 10000;
  // NOTE TODOL current implementation requires explicit touch-create before detailed balance
  map<array<double, Specs>> particles;
  {
    auto& x = particles[species::electron];
    auto& y = particles[species::ion];
    auto& z = particles[species::photon];
  }

  auto rwld_opt = aio::reduced_world(num_procs, mpi::world);
  if ( rwld_opt ) {
    auto cart_opt = aio::make_cart( {num_ens}, {false}, *rwld_opt );
    auto ens_opt = dye::create_ensemble<1>(cart_opt);

    aio::gauss_real<double> load_gen( 1000.0, 500.0 );
    load_gen.seed( aio::now() + rwld_opt->rank() );
    unsigned int target_load = 0;

    int N = num_cycles;
    while ( N-- ) {
      if ( ens_opt ) {
        for ( auto&[sp, ptcs] : particles ) {
            auto x = load_gen();
            x = std::max<double>(x, 0.0);
            x = std::min<double>(x, 100000.0);
            ptcs.resize( static_cast<load_t>(x) );
        }
      }
      dye::dynamic_load_balance(particles, ens_opt, cart_opt, target_load );
    }

  }
  mpi::world.barrier();

}
