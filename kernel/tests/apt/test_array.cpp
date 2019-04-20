#include "apt/array.hpp"
#include "catch2/catch.hpp"

using namespace apt;

SCENARIO("array deepcopy constructor", "[apt]" ) {
  array<int,3> x{ 3, 4, 7 };
  array<int,3> arr(x);
  REQUIRE_FALSE( x._data == arr._data );
  x[0] = 213;
  x[1] = 213;
  x[2] = 213;
  REQUIRE( arr[0] == 3 );
  REQUIRE( arr[1] == 4 );
  REQUIRE( arr[2] == 7 );
}

SCENARIO("array deepcopy assign", "[apt]" ) {
  array<int,3> x{ 3, 4, 7 };
  array<int,3> arr{ 2, 6, 8 };
  arr = x;
  REQUIRE_FALSE( x._data == arr._data );
  x[0] = 213;
  x[1] = 213;
  x[2] = 213;
  REQUIRE( arr[0] == 3 );
  REQUIRE( arr[1] == 4 );
  REQUIRE( arr[2] == 7 );
}

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

SCENARIO("nested array", "[apt]" ) {
  array< array<int,2>, 3 > arrs;
  for ( int j = 0; j < 3; ++j )
    for ( int i = 0; i < 2; ++i )
      REQUIRE( arrs[j][i] == 0 );
}

