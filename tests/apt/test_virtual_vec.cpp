#include "apt/virtual_vec.hpp"
#include "catch2/catch.hpp"

using vVec = apt::vVec<double,3>;

SCENARIO( "vVec constructors", "[apt]" ) {
  SECTION("value constructor and reference sematics") {
    double x = 1.2, y = 2.3;
    vVec v ( x, y, y );
    REQUIRE( std::get<0>(v) == 1.2 );
    REQUIRE( std::get<1>(v) == 2.3 );
    REQUIRE( std::get<2>(v) == 2.3 );

    std::get<0>(v) = 12.0;
    std::get<2>(v) = 23.0;

    REQUIRE( x == 12.0 );
    REQUIRE( y == 23.0 );
    REQUIRE( std::get<0>(v) == 12.0 );
    REQUIRE( std::get<1>(v) == 23.0 );
    REQUIRE( std::get<2>(v) == 23.0 );
  }

  SECTION("move constructors") {
    double x = 1.2, y = 2.3, z = 3.4;
    vVec v0 ( x, y, z );
    vVec v_mv = std::move(v0);
    REQUIRE( std::get<0>(v_mv) == 1.2 );
    REQUIRE( std::get<1>(v_mv) == 2.3 );
    REQUIRE( std::get<2>(v_mv) == 3.4 );
  }
}
