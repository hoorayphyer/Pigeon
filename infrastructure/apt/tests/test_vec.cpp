#include "apt/vec.hpp"
#include "apt/print.hpp"
#include "catch2/catch.hpp"

using Vec = apt::Vec<double,3>;

SCENARIO( "Vec constructors", "[apt][vec]" ) {
  SECTION("default constructor") {
    Vec v;
    REQUIRE( v[0] == 0.0 );
    REQUIRE( v[1] == 0.0 );
    REQUIRE( v[2] == 0.0 );
  }

  SECTION("value constructor") {
    Vec v ( 1.2, 2.3, 3.4 );
    REQUIRE( v[0] == 1.2 );
    REQUIRE( v[1] == 2.3 );
    REQUIRE( v[2] == 3.4 );
  }

  SECTION("copy/move constructors") {
    Vec v0 ( 1.2, 2.3, 3.4 );
    Vec v_cp = v0;
    REQUIRE( v_cp[0] == 1.2 );
    REQUIRE( v_cp[1] == 2.3 );
    REQUIRE( v_cp[2] == 3.4 );

    Vec v_mv = std::move(v0);
    REQUIRE( v_mv[0] == 1.2 );
    REQUIRE( v_mv[1] == 2.3 );
    REQUIRE( v_mv[2] == 3.4 );
  }
}

SCENARIO( "Vec lvalue ref accessor", "[apt][vec]") {
  double x = 1.2;
  Vec v ( x, x, x );
  v[1] = 3.4;

  REQUIRE( x == 1.2 );
  REQUIRE( v[0] == 1.2 );
  REQUIRE( v[1] == 3.4 );
  REQUIRE( v[2] == 1.2 );
}
