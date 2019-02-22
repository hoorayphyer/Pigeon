#include "apt/vec.hpp"
#include "catch2/catch.hpp"

using Vec = apt::Vec<double,3>;

SCENARIO( "Vec constructors", "[apt]" ) {
  SECTION("default constructor") {
    Vec v;
    REQUIRE( std::get<0>(v) == 0.0 );
    REQUIRE( std::get<1>(v) == 0.0 );
    REQUIRE( std::get<2>(v) == 0.0 );
  }

  SECTION("value constructor") {
    Vec v ( 1.2, 2.3, 3.4 );
    REQUIRE( std::get<0>(v) == 1.2 );
    REQUIRE( std::get<1>(v) == 2.3 );
    REQUIRE( std::get<2>(v) == 3.4 );
  }

  SECTION("copy/move constructors") {
    Vec v0 ( 1.2, 2.3, 3.4 );
    Vec v_cp = v0;
    REQUIRE( std::get<0>(v_cp) == 1.2 );
    REQUIRE( std::get<1>(v_cp) == 2.3 );
    REQUIRE( std::get<2>(v_cp) == 3.4 );

    Vec v_mv = std::move(v0);
    REQUIRE( std::get<0>(v_mv) == 1.2 );
    REQUIRE( std::get<1>(v_mv) == 2.3 );
    REQUIRE( std::get<2>(v_mv) == 3.4 );
  }
}

SCENARIO( "Vec lvalue ref accessor", "[apt]") {
  double x = 1.2;
  Vec v ( x, x, x );
  std::get<1>(v) = 3.4;

  REQUIRE( x == 1.2 );
  REQUIRE( std::get<0>(v) == 1.2 );
  REQUIRE( std::get<1>(v) == 3.4 );
  REQUIRE( std::get<2>(v) == 1.2 );
}
