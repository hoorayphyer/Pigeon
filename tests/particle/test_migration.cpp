#include "all_in_one.hpp"
#include "particle/migration.cpp"
#include "particle/particle.hpp"
#include "dye/ensemble.hpp"

using namespace particle;
using cPtc_t = cParticle<double,3,unsigned long long>;
using Ptc_t = Particle<double,3,unsigned long long>;
using Vec_t = apt::Vec<double,3>;

SCENARIO("Test sendrecv cparticle", "[particle][mpi][.]") {
  if ( mpi::world.size() > 1  ) {
    aio::unif_real<double> dist;
    std::vector<Ptc_t> ptc_arr;
    const int nptcs = 1000;
    ptc_arr.reserve(nptcs);
    for ( int i = 0; i < nptcs; ++i )
      ptc_arr.emplace_back( Vec_t(dist(),dist(),dist()), Vec_t{dist(),dist(),dist()}, flag::secondary, species::electron );

    std::vector<cPtc_t> data;
    data.reserve(nptcs);
    data.resize(nptcs);

    if ( mpi::world.rank() == 0 ) {
      for ( int i = 0; i < nptcs; ++i )
        data[i] = ptc_arr[i];

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

}

SCENARIO("Test lcr_sort", "[particle][.]") {
  if ( mpi::world.rank() == 0 ) {
    auto test =
      []( const std::vector<int>& lcr_vals ) {

        auto f_lcr = []( const auto& ptc ) { return impl::lcr(ptc.q()[0], 0.5, 1.5); };
        std::vector<cPtc_t> ptcs;

        apt::array<int,3> count{};

        for ( int i = 0; i < lcr_vals.size(); ++i ) {
          ++count[ lcr_vals[i] % 3 ];
          Ptc_t p;
          p.q()[0] = static_cast<double>( lcr_vals[i] % 3 );
          ptcs.emplace_back(std::move(p));
        }

        auto begs = impl::lcr_sort( ptcs, f_lcr );

        int begL_exp = count[1];
        int begR_exp = count[0] + count[1];
        int begE_exp = count[0] + count[1] + count[2];

        REQUIRE( begs[0] == begL_exp );
        REQUIRE( begs[1] == begR_exp );
        REQUIRE( begs[2] == begE_exp );

        int n = 0;
        for ( ; n < begs[0]; ++n ) REQUIRE( ptcs[n].q()[0] == 1.0 );
        for ( ; n < begs[1]; ++n ) REQUIRE( ptcs[n].q()[0] == 0.0 );
        for ( ; n < begs[2]; ++n ) REQUIRE( ptcs[n].q()[0] == 2.0 );

      };

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

    // to uniformly sample the peripheries of the bulk with high hit rate, we sample a square torus. Points are generated in this range but those that fall within the inner square are dumped. Since no action is performed after particle migrates, it is OK to have particle position values that jump across whole node.
    constexpr double outer = 5.0;
    constexpr double inner = 1.0;

    apt::array< apt::pair<double>, DGrid> borders;
    for ( auto& b : borders ) b = { - inner, inner };

    aio::unif_real<double> dist( - outer, outer );
    dist.seed( aio::now() + mpi::world.rank() );

    int nptcs = 10000;
    constexpr int N = 3*3; // 3^DGrid, including self
    constexpr int CENTER = ( N - 1 ) / 2;
    apt::array<int, N> sendcount{};
    while( nptcs ) {
      apt::array<double, DGrid> smp;
      bool in_bulk = true;
      for ( int i = 0; i < DGrid; ++i ) {
        smp[i] = dist();
        in_bulk = ( in_bulk && std::abs(smp[i]) < inner );
      }
      if (in_bulk) continue;

      Ptc_t ptc;
      for ( int i = 0; i < DGrid; ++i ) {
        ptc.q()[i] = smp[i];
        ptc.p()[i] = ens_opt -> cart_coords[i]; // use p for encoding
      }
      CAPTURE( ptc.q(), borders );
      REQUIRE( is_migrate(ptc.q(), borders) );

      mgr_buf.push_back( std::move(ptc) );
      ++sendcount[ impl::lcr(smp[0], -inner, inner) + 3 * impl::lcr(smp[1], -inner, inner) ];
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

    migrate( mgr_buf, inter, borders, 0u );

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
