#include "apt/array.hpp"
#include "catch2/catch.hpp"

using namespace apt;

SCENARIO("array comparison", "[apt]") {
  array<int,3> x{ 3, 4, 7 };
  array<int,3> y{ 3, 7, 1 };
  auto z (x);

  REQUIRE( x == x );
  REQUIRE_FALSE( x == y );
  REQUIRE( x == z );
}
