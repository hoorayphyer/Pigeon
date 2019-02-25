#include "apt/virtual_vec.hpp"
#include "catch2/catch.hpp"

using vVec = apt::vVec<double,3>;

SCENARIO( "vVec constructors", "[apt][vec]" ) {
  SECTION("value constructor and reference sematics") {
    double x = 1.2, y = 2.3;
    vVec v ( x, y, y );
    REQUIRE( v[0] == 1.2 );
    REQUIRE( v[1] == 2.3 );
    REQUIRE( v[2] == 2.3 );

    v[0] = 12.0;
    v[2] = 23.0;

    REQUIRE( x == 12.0 );
    REQUIRE( y == 23.0 );
    REQUIRE( v[0] == 12.0 );
    REQUIRE( v[1] == 23.0 );
    REQUIRE( v[2] == 23.0 );
  }

  SECTION("constructor from std::array") {
    std::array<double,3> arr { 1, 2, 3 };
    vVec v (arr);

    REQUIRE( v[0] == 1 );
    REQUIRE( v[1] == 2 );
    REQUIRE( v[2] == 3 );
  }

  SECTION("move constructors") {
    double x = 1.2, y = 2.3, z = 3.4;
    vVec v0 ( x, y, z );
    vVec v_mv = std::move(v0);
    REQUIRE( v_mv[0] == 1.2 );
    REQUIRE( v_mv[1] == 2.3 );
    REQUIRE( v_mv[2] == 3.4 );
  }
}
