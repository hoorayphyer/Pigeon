#include "aperture/ensemble.cpp"
#include "parallel/mpi_datatype.cpp"
#include "parallel/mpi++.cpp"
#include "parallel/mpi_communication.cpp"
#include "catch2/catch.hpp"

using namespace aperture;
using namespace mpi;

// TODO do a larger test
// TODO add more tests
SCENARIO("Create Ensemble", "[parallel][aperture]") {
  constexpr int DGrid = 1;
  if ( world.size() == 2 ) {
    SECTION("one proc as the primary, together with the other it forms an ensemble") {
      auto cart_opt = create_primary_comm({1}, {false});
      if ( world.rank() == 0 ) {
        REQUIRE( cart_opt );
      } else {
        REQUIRE_FALSE( cart_opt );
      }
      auto intra_opt = world.split(Group({0,1}));
      REQUIRE(intra_opt);

      auto ens_opt = create_ensemble<DGrid>( cart_opt, intra_opt );
      REQUIRE( ens_opt );
      for ( int i = 0; i < DGrid; ++i )
        for ( int b = 0; b < 2; ++b )
          REQUIRE_FALSE( ens_opt->inter[i][b] );
      const auto& ens = *ens_opt;
      REQUIRE( 0 == ens.chief );
      REQUIRE( 0 == ens.label );
      REQUIRE( 0 == ens.chief_cart_rank );

      REQUIRE( apt::array<int,1>{0} == ens.cart_coords );
      for ( int IGrid = 0; IGrid < DGrid; ++ IGrid ) {
        for ( int i = 0; i < 2; ++i ) {
          REQUIRE_FALSE( ens.neigh_cart_ranks[IGrid][i] );
          REQUIRE( ens.is_at_boundary()[IGrid][i] );
        }
      }
    }
  }
}

SCENARIO("Link Neighbors", "[parallel]") {
  constexpr int DGrid = 1;
  if ( world.size() == 2 ) {
    auto cart_opt = create_primary_comm({2}, {true});
    REQUIRE( cart_opt );
    auto intra_opt = world.split(Group({world.rank()}));
    REQUIRE( intra_opt );
    auto ens_opt = create_ensemble<DGrid>( cart_opt, intra_opt );
    REQUIRE( ens_opt );

    // TODO
    // auto neighbors = link_neighbors( cart_opt, intra_opt, locale_opt )[0];

    // const auto& inter = ens_opt->inter[0];
    // for ( int i = 0; i < 2; ++i ) {
    //   REQUIRE( inter[i] );
    //   const auto& itc = *(inter[i]);
    //   REQUIRE( itc.remote_size() == 1 );
    //   // REQUIRE( itc.remote_group() ==  ); // TODO group comparison
    // }
  }
}
