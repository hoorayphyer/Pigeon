#include "parallel/locale.cpp"
#include "parallel/mpi++.cpp"
#include "parallel/mpi_communication.cpp"
#include "catch2/catch.hpp"

using namespace parallel;
using namespace mpi;

// TODO do a larger test
SCENARIO("Create Locale", "[parallel]") {
  constexpr int DGrid = 1;
  if ( world.size() == 2 ) {
    SECTION("one proc as the primary, together with the other it forms an ensemble") {
      auto cart_opt = create_primary_comm({1}, {false});
      if ( world.rank() == 0 ) {
        REQUIRE( cart_opt );
      } else {
        REQUIRE_FALSE( cart_opt );
      }
      auto ens_opt = world.split(Group({0,1}));
      REQUIRE(ens_opt);

      auto locale = create_locale<DGrid>( cart_opt, *ens_opt );
      REQUIRE( 0 == locale.chief );
      REQUIRE( 0 == locale.label );
      REQUIRE( 0 == locale.chief_cart_rank );

      REQUIRE( std::array<int,1>{0} == locale.cart_coords );
      for ( int IGrid = 0; IGrid < DGrid; ++ IGrid )
        for ( int i = 0; i < 2; ++i ) {
          REQUIRE_FALSE( locale.neighbors[IGrid][i] );
          REQUIRE( locale.is_at_boundary[IGrid][i] );
        }
    }

  }
}

SCENARIO("Link Neighbors", "[parallel]") {
  
}
