#include "apt/pair.hpp"
#include "catch2/catch.hpp"

using namespace apt;

SCENARIO("pair", "[apt]") {
  pair<int> x{ 6, 8 };
  REQUIRE( x[0] == 6 );
  REQUIRE( x[LFT] == 6 );
  REQUIRE( x[1] == 8 );
  REQUIRE( x[RGT] == 8 );
  x[0]++;
  REQUIRE( x[0] == 7 );
}
