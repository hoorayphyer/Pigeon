#include "testfw/testfw.hpp"
#include "pic.hpp"
#include "particle/migration_impl.hpp"
#include "particle/particle.hpp"
#include "dye/ensemble_impl.hpp"

using namespace particle;
using Real = pic::real_t;
using Ptc_t = Particle<Real,Specs>;
using Vec_t = apt::Vec<Real,Specs<Real>::Dim>;

SCENARIO("Test swapping particles", "[particle]") {
  if ( mpi::world.rank() == 0 ) {
    Ptc_t a( {3.5,3.5,3.5 }, {5.9, 5.9, 5.9}, flag::secondary, birthplace{3} );
    migrInt<2>(6).imprint(a);

    Ptc_t b( {-4.7, -4.7, -4.7, }, {-2.8, -2.8, -2.8}, birthplace{98} );
    migrInt<2>(1).imprint(b);
    b.reset(flag::exist);

    a.swap(b);

    {
      REQUIRE( a.q()[0] == static_cast<Real>(-4.7) );
      REQUIRE( a.q()[1] == static_cast<Real>(-4.7) );
      REQUIRE( a.q()[2] == static_cast<Real>(-4.7) );

      REQUIRE( a.p()[0] == static_cast<Real>(-2.8) );
      REQUIRE( a.p()[1] == static_cast<Real>(-2.8) );
      REQUIRE( a.p()[2] == static_cast<Real>(-2.8) );

      REQUIRE_FALSE( a.is(flag::exist) );
      REQUIRE_FALSE( a.is(flag::secondary) );
      REQUIRE( a.get<birthplace>() == 98 );
      REQUIRE( migrInt<2>(a) == 1 );
    }

    {
      REQUIRE( b.q()[0] == static_cast<Real>(3.5) );
      REQUIRE( b.q()[1] == static_cast<Real>(3.5) );
      REQUIRE( b.q()[2] == static_cast<Real>(3.5) );

      REQUIRE( b.p()[0] == static_cast<Real>(5.9) );
      REQUIRE( b.p()[1] == static_cast<Real>(5.9) );
      REQUIRE( b.p()[2] == static_cast<Real>(5.9) );

      REQUIRE( b.is(flag::exist) );
      REQUIRE( b.is(flag::secondary) );
      REQUIRE( b.get<birthplace>() == 3 );
      REQUIRE( migrInt<2>(b) == 6 );
    }

  }
}

SCENARIO("Test sending particles", "[particle][mpi]") {
  SECTION("Intra communicator") {
    if ( mpi::world.size() >= 2 ) {
      int myrank = mpi::world.rank();
      if ( 0 == myrank ) {
        Ptc_t ptc( {1.1,2.2,3.3}, {4.4,5.5,6.6}, birthplace{18} );
        mpi::world.send( 1, 22, &ptc );
      } else if ( 1 == myrank) {
        Ptc_t ptc;
        mpi::world.recv( 0, 22, &ptc, 1 );

        REQUIRE( ptc.q()[0] == static_cast<Real>(1.1) );
        REQUIRE( ptc.q()[1] == static_cast<Real>(2.2) );
        REQUIRE( ptc.q()[2] == static_cast<Real>(3.3) );

        REQUIRE( ptc.p()[0] == static_cast<Real>(4.4) );
        REQUIRE( ptc.p()[1] == static_cast<Real>(5.5) );
        REQUIRE( ptc.p()[2] == static_cast<Real>(6.6) );

        REQUIRE( ptc.is(flag::exist) );
        REQUIRE( ptc.get<birthplace>() == 18 );
      }
    }
  }

  SECTION("Inter communicator") {
    if ( mpi::world.size() >= 2 && mpi::world.rank() < 2 ) {
      int myrank = mpi::world.rank();
      mpi::InterComm inter( mpi::self, 0, {mpi::world}, 1-myrank, 147 );
      if ( 0 == myrank ) {
        Ptc_t ptc( {1.1,2.2,3.3}, {4.4,5.5,6.6}, birthplace{18} );
        inter.send( 0, 22, &ptc );
      } else if ( 1 == myrank) {
        Ptc_t ptc;
        inter.recv( 0, 22, &ptc, 1 );

        REQUIRE( ptc.q()[0] == static_cast<Real>(1.1) );
        REQUIRE( ptc.q()[1] == static_cast<Real>(2.2) );
        REQUIRE( ptc.q()[2] == static_cast<Real>(3.3) );

        REQUIRE( ptc.p()[0] == static_cast<Real>(4.4) );
        REQUIRE( ptc.p()[1] == static_cast<Real>(5.5) );
        REQUIRE( ptc.p()[2] == static_cast<Real>(6.6) );

        REQUIRE( ptc.is(flag::exist) );
        REQUIRE( ptc.get<birthplace>() == 18 );
      }
    }
  }
}

TEMPLATE_TEST_CASE("Test setting destination", "[particle]"
                   , (std::integral_constant<int,2>)
                   , (std::integral_constant<int,3>)
                   ) {
  constexpr int DGrid = TestType::value;
  Ptc_t ptc;
  REQUIRE( migrInt<DGrid>(ptc) == ( apt::pow3(DGrid) - 1 ) / 2 );

  for ( int i = 0; i < apt::pow3(DGrid); ++i ) {
    migrInt<DGrid> mint(i);
    mint.imprint(ptc);
    REQUIRE( migrInt<DGrid>(ptc) == i );
  }
}

SCENARIO("Test sendrecv particle", "[particle][mpi]") {
  if ( mpi::world.size() > 1 && mpi::world.rank() < 2  ) {
    aio::unif_real<double> dist;
    std::vector<Ptc_t> ptc_arr;
    const int nptcs = 1000;
    ptc_arr.reserve(nptcs);
    for ( int i = 0; i < nptcs; ++i ) {
      ptc_arr.emplace_back( Vec_t(dist(),dist(),dist()), Vec_t{dist(),dist(),dist()}, flag::secondary, species::electron );
      migrInt<3>(14).imprint(ptc_arr.back());
    }

    std::vector<Ptc_t> data;
    data.reserve(nptcs);
    data.resize(nptcs);
    for ( auto& ptc : data ) ptc.set(flag::exist);

    if ( mpi::world.rank() == 0 ) {
      for ( int i = 0; i < nptcs; ++i ) {
        data[i] = ptc_arr[i];
        migrInt<3>(14).imprint(data[i]);
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
    }
  }
  mpi::world.barrier();
}

SCENARIO("Test lcr_sort", "[particle]") {
  auto test =
    []( const std::vector<int>& lcr_vals ) {
      if ( mpi::world.rank() == 0 ) {

        auto lcr = []( const auto& ptc ) -> int { return migrInt<1>(ptc) % 3; };
        std::vector<Ptc_t> ptcs;

        apt::array<int,3> count{};

        for ( int i = 0; i < lcr_vals.size(); ++i ) {
          ++count[ lcr_vals[i] % 3 ];
          ptcs.emplace_back(Ptc_t());
          migrInt<1>(lcr_vals[i] % 3).imprint(ptcs.back());
        }

        auto begs = impl::lcr_sort( ptcs, lcr );

        int begL_exp = count[1];
        int begR_exp = count[0] + count[1];
        int begE_exp = count[0] + count[1] + count[2];

        REQUIRE( begs[0] == begL_exp );
        REQUIRE( begs[1] == begR_exp );
        REQUIRE( begs[2] == begE_exp );

        int n = 0;
        for ( ; n < begs[0]; ++n ) REQUIRE( migrInt<1>(ptcs[n]) == apt::C );
        for ( ; n < begs[1]; ++n ) REQUIRE( migrInt<1>(ptcs[n]) == apt::L );
        for ( ; n < begs[2]; ++n ) REQUIRE( migrInt<1>(ptcs[n]) == apt::R );

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

TEMPLATE_TEST_CASE("Testing particle migration with trivial ensemble", "[field][mpi]"
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

  auto cart_opt = aio::make_cart( TestType::get(), mpi::world );
  auto ens_opt = dye::create_ensemble<DGrid>(cart_opt);
  if ( ens_opt ) {
    // first every node initialize some particles out of borders. Then they share this information with neighbors for test purpose.
    std::vector<Ptc_t> mgr_buf;

    aio::unif_int<int> unif( 0, apt::pow3(DGrid) - 1 );
    unif.seed( aio::now() + mpi::world.rank() );

    int nptcs = 10000;
    constexpr int N = apt::pow3(DGrid);
    constexpr int CENTER = ( N - 1 ) / 2;
    apt::array<int, N> sendcount{};
    while( nptcs ) {
      migrInt<DGrid> mig_dir = unif();
      if (mig_dir == CENTER) continue;
      ++sendcount[mig_dir];

      Ptc_t ptc;
      ptc.set(flag::exist);
      for ( int i = 0; i < DGrid; ++i ) {
        ptc.p()[i] = ens_opt -> cart_coords[i]; // use p for encoding
      }

      mig_dir.imprint(ptc);
      mgr_buf.push_back( std::move(ptc) );
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

    std::vector<Ptc_t> ptcs;
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
          {
            auto x = load_gen();
            x = std::max<double>(x, 0.0);
            x = std::min<double>(x, 100000.0);
            ptcs.resize(static_cast<std::size_t>(x));
          }
          for ( auto& ptc : ptcs ) {
            ptc.set(flag::exist);
            migrInt<DGrid> mi;
            do mi = mig_dir_gen();
            while ( mi == (apt::pow3(DGrid) - 1) / 2 );
            mi.imprint(ptc);
          }
        }
        migrate( ptcs, ens_opt->inter, 0 ); // TODO shift
      }
    }

  }
  mpi::world.barrier();
}
