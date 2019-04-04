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

SCENARIO("range loop", "[apt]") {
  array<int,6> x{ 3, 4, 7, 5, 6, 2 };
  auto z (x);

  int i = 0;
  for ( const auto& elm : x ) {
    REQUIRE( elm == z[i++] );
  }
  REQUIRE( i == x.size() );
}
