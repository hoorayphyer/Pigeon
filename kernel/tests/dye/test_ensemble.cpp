#include "all_in_one.hpp"
#include "dye/ensemble.hpp"

TEMPLATE_TEST_CASE( "Create Ensemble","[dye][mpi]"
                   // NOTE Notation: XxYxZ is the cartesian partition. The cartesian topology is periodic in all directions
                   , (aio::IndexType<-1,-1>)
                   , (aio::IndexType<-2,-1>)
                   , (aio::IndexType<-1,-2>)
                   , (aio::IndexType<-2,-2>)
                   , (aio::IndexType<-4,-4>)
                   , (aio::IndexType<-8,-8>)
                   ) {
  constexpr int DGrid = 2;
  std::vector<int> cart_dims;
  std::vector<bool> periodic;
  for ( auto i : TestType::get() ) {
    cart_dims.push_back( i > 0 ? i : -i );
    periodic.push_back( true );
  }

  WHEN("creating trivial ensemble") {
    auto cart_opt = aio::make_cart( cart_dims, periodic, mpi::world );
    auto ens_opt = dye::create_ensemble<DGrid>(cart_opt);
    REQUIRE( bool(cart_opt) == bool(ens_opt) );
    if ( ens_opt ) {
      const auto& cart = *cart_opt;
      const auto& ens = *ens_opt;
      REQUIRE(ens.intra.size() == 1);
      for ( int i = 0; i < DGrid; ++i ) {
        auto[src,dest] = cart.shift(i);
        if ( 1 == cart_dims[i] && periodic[i] ) {
          REQUIRE_FALSE( ens.inter[i][LFT] );
          REQUIRE_FALSE( ens.inter[i][RGT] );
        } else {
          REQUIRE( bool(src) == bool(ens.inter[i][LFT]) );
          REQUIRE( bool(dest) == bool(ens.inter[i][RGT]) );
        }
      }

      REQUIRE(ens.chief == 0);
      REQUIRE(ens.chief_cart_rank == cart.rank());
      auto[ c, d, p ] = cart.coords_dims_periodic();
      for ( int i = 0; i < DGrid; ++i ) {
        REQUIRE( ens.cart_coords[i] == c[i] );
        REQUIRE( ens.cart_dims[i] == d[i] );
        REQUIRE( ens.is_periodic[i] == p[i] );
      }
    }
  }
}
