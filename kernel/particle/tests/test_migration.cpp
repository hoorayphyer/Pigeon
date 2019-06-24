#include "testfw/testfw.hpp"
#include "pic.hpp"
#include "particle/migration_impl.hpp"
#include "particle/particle.hpp"
#include "dye/ensemble_impl.hpp"

using namespace particle;
using cPtc_t = cParticle<double,Specs>;
using Ptc_t = Particle<double,Specs>;
using Vec_t = apt::Vec<double,Specs<double>::Dim>;

SCENARIO("Test sendrecv cparticle", "[particle][mpi][.]") {
  if ( mpi::world.size() > 1 && mpi::world.rank() < 2  ) {
    aio::unif_real<double> dist;
    std::vector<Ptc_t> ptc_arr;
    const int nptcs = 1000;
    ptc_arr.reserve(nptcs);
    for ( int i = 0; i < nptcs; ++i )
      ptc_arr.emplace_back( Vec_t(dist(),dist(),dist()), Vec_t{dist(),dist(),dist()}, flag::secondary, species::electron );

    std::vector<cPtc_t> data;
    data.reserve(nptcs);
    data.resize(nptcs);
    for ( auto& ptc : data ) ptc.set(flag::exist);

    if ( mpi::world.rank() == 0 ) {
      for ( int i = 0; i < nptcs; ++i ) {
        data[i] = ptc_arr[i];
        data[i].extra() = 14;
      }

      mpi::world.send(1,147,data.data(),data.size());
    } else if ( mpi::world.rank() == 1 ) {
      mpi::world.recv(0,147,data.data(),data.size());
    }

    for ( int i = 0; i < nptcs; ++i ) {
      const auto& cptc = data[i];
      const auto& ptc = ptc_arr[i];
      for ( int n = 0; n < 3; ++n ) {
        REQUIRE( cptc.q()[n] == ptc.q()[n] );
        REQUIRE( cptc.p()[n] == ptc.p()[n] );
      }
      REQUIRE( cptc.state() == ptc.state() );
      REQUIRE( cptc.extra() == 14 );
    }
  }
  mpi::world.barrier();
}

SCENARIO("Test lcr_sort", "[particle][.]") {
  auto test =
    []( const std::vector<int>& lcr_vals ) {
      if ( mpi::world.rank() == 0 ) {

        auto lcr = []( const auto& ptc ) { return ptc.extra() % 3; };
        std::vector<cPtc_t> ptcs;

        apt::array<int,3> count{};

        for ( int i = 0; i < lcr_vals.size(); ++i ) {
          ++count[ lcr_vals[i] % 3 ];
          ptcs.emplace_back(Ptc_t());
          ptcs.back().extra() = lcr_vals[i] % 3;
        }

        auto begs = impl::lcr_sort( ptcs, lcr );

        int begL_exp = count[1];
        int begR_exp = count[0] + count[1];
        int begE_exp = count[0] + count[1] + count[2];

        REQUIRE( begs[0] == begL_exp );
        REQUIRE( begs[1] == begR_exp );
        REQUIRE( begs[2] == begE_exp );

        int n = 0;
        for ( ; n < begs[0]; ++n ) REQUIRE( ptcs[n].extra() == 1 );
        for ( ; n < begs[1]; ++n ) REQUIRE( ptcs[n].extra() == 0 );
        for ( ; n < begs[2]; ++n ) REQUIRE( ptcs[n].extra() == 2 );

      }
      mpi::world.barrier();
    };

  // NOTE seems that all multiple each SECTION will rerun the whole test, so if there is mpi barrier needed, all processes have to enter the SECTION to avoid hanging behavior
  SECTION("empty buffer") {
    test(std::vector<int>{});
  }

  SECTION("uniform values") {
    for ( int v = 0; v < 3; ++v ) {
      const int nptcs = 100;
      std::vector<int> lcr_vals(nptcs);
      for ( auto& x : lcr_vals ) x = v;
      test(lcr_vals);
    }
  }

  SECTION("random test") {
    aio::unif_int<int> unif(0, 2);
    int N = 10;
    const int nptcs = 1000;
    while ( N-- ) {
      std::vector<int> lcr_vals(nptcs);
      for ( auto& x : lcr_vals ) x = unif();
      test(lcr_vals);
    }
  }
}

TEMPLATE_TEST_CASE("Testing particle migration with trivial ensemble", "[field][mpi][.]"
                   // NOTE Notation: XxYxZ is the cartesian partition. The cartesian topology is periodic in all directions
                   , (aio::IndexType<1,1>)
                   , (aio::IndexType<2,1>)
                   , (aio::IndexType<1,2>)
                   , (aio::IndexType<2,2>)
                   // , (aio::IndexType<-1,-1>)
                   // , (aio::IndexType<-2,-1>)
                   // , (aio::IndexType<-1,-2>)
                   // , (aio::IndexType<-2,-2>)
                   // , (aio::IndexType<-4,-4>)
                   // , (aio::IndexType<-8,-8>)
                   ) {
  constexpr int DGrid = 2;
  std::vector<int> cart_dims;
  std::vector<bool> periodic;
  for ( auto i : TestType::get() ) {
    cart_dims.push_back( i > 0 ? i : -i );
    periodic.push_back( i < 0 );
  }

  auto cart_opt = aio::make_cart( cart_dims, periodic, mpi::world );
  auto ens_opt = dye::create_ensemble<DGrid>(cart_opt);
  if ( ens_opt ) {
    // first every node initialize some particles out of borders. Then they share this information with neighbors for test purpose.
    std::vector<cPtc_t> mgr_buf;

    aio::unif_int<int> unif( 0, apt::pow3(DGrid) - 1 );
    unif.seed( aio::now() + mpi::world.rank() );

    int nptcs = 10000;
    constexpr int N = apt::pow3(DGrid);
    constexpr int CENTER = ( N - 1 ) / 2;
    apt::array<int, N> sendcount{};
    while( nptcs ) {
      int mig_dir = unif();
      if (mig_dir == CENTER) continue;
      ++sendcount[mig_dir];

      Ptc_t ptc;
      for ( int i = 0; i < DGrid; ++i ) {
        ptc.p()[i] = ens_opt -> cart_coords[i]; // use p for encoding
      }

      mgr_buf.push_back( std::move(ptc) );
      mgr_buf.back().extra() = mig_dir;
      -- nptcs;
    }

    apt::array<int, N> recvcount{};
    std::vector<mpi::Request> reqs;
    const auto[ my_coords, cart_dims, periodic ] = cart_opt -> coords_dims_periodic();
    for( int n = 0; n < N; ++n ) {
      if ( CENTER == n ) continue;
      auto neigh_coords = my_coords;
      neigh_coords[1] += n / 3 - 1;
      neigh_coords[0] += n % 3 - 1;

      std::optional<int> neigh_rank = cart_opt -> coords2rank(neigh_coords);
      if ( neigh_rank ) {
        reqs.push_back( cart_opt -> Isend(*neigh_rank, 147, sendcount.begin() + n, 1) );
        reqs.push_back( cart_opt -> Irecv(*neigh_rank, 147, recvcount.begin() + n, 1) );
      }
    }
    mpi::waitall(reqs);

    const auto& inter = ens_opt -> inter;

    migrate( mgr_buf, inter, 0u );

    apt::array<int,N> recv_actual{};
    for ( const auto& ptc : mgr_buf ) {
      auto neigh_rel = my_coords;
      for ( int i = 0; i < DGrid; ++i)
        neigh_rel[i] = std::lround(ptc.p()[i]) - my_coords[i] + 1; // neigh_rel = 0, 1, 2

      ++recv_actual[ neigh_rel[0] + 3 * neigh_rel[1] ];
    }

    CAPTURE( my_coords, sendcount, recvcount, recv_actual );
    for( int i = 0; i < N; ++i ) {
      REQUIRE( recv_actual[i] == recvcount[i] );
    }

  }
  mpi::world.barrier();
}

SCENARIO("Stress Test" , "[particle][mpi]") {
  constexpr int DGrid = 2;
  const apt::array<int,DGrid> partition { 3, 3 };
  const int num_procs = mpi::world.size();
  const int num_proc_map = 1000;
  const int num_cycles_each_proc_map = 10000;

  auto rwld_opt = aio::reduced_world(num_procs, mpi::world);
  if ( rwld_opt ) {
    auto cart_opt = aio::make_cart( partition, *rwld_opt );
    int M = num_proc_map;
    int num_ens = 1;
    for (auto x : partition) num_ens *= x;
    // TODO seeding
    aio::unif_int<int> ui( 0, num_ens - 1 );
    aio::gauss_real<double> load_gen( 1000.0, 1000.0 );
    // load_gen.seed( aio::now() + rwld_opt->rank() );
    load_gen.seed( rwld_opt->rank() );
    aio::unif_int<int> mig_dir_gen( 0, apt::pow3(DGrid) - 1 );

    std::vector<cPtc_t> ptcs;
    while ( M-- ) {
      int color = 0;
      if ( rwld_opt->rank() < num_ens ) color = rwld_opt->rank();
      else color = ui();
      std::optional<dye::Ensemble<DGrid>> ens_opt;
      {
        auto intra = rwld_opt->split(color);
        ens_opt = dye::create_ensemble<DGrid>(cart_opt, intra);
      }
      int N = num_cycles_each_proc_map;
      while ( N-- ) {
        {
          auto x = load_gen();
          x = std::max<double>(x, 0.0);
          x = std::min<double>(x, 100000.0);
          ptcs.resize(static_cast<std::size_t>(x));
          for ( auto& ptc : ptcs ) {
            do {
              ptc.set(flag::exist);
              ptc.extra() = mig_dir_gen();
            } while ( ptc.extra() == (apt::pow3(DGrid) - 1) / 2 );
          }
        }
        migrate( ptcs, ens_opt->inter, 0 ); // TODO shift
      }
    }

  }
  mpi::world.barrier();
}
