#include "apt/block.hpp"
#include "apt/print.hpp"
#include "catch2/catch.hpp"

using namespace apt;

SCENARIO("Block", "[apt]") {
  int N[3] = { 13, 15, 17 };
  Block block{ Index<3>{ N[0], N[1], N[2] }};
  auto itr = block.begin();
  auto end = block.end();
  REQUIRE( end == Index<3>{ N[0], N[1]-1, N[2]-1 } );

  for ( int k = 0; k < N[2]; ++k )
    for ( int j = 0; j < N[1]; ++j )
      for ( int i = 0; i < N[0]; ++i ) {
        CAPTURE(*itr, i, j, k);
        REQUIRE( *(itr++) == Index<3>{ i, j, k } );
      }
}
