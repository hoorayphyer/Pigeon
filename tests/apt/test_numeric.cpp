#include "apt/numeric.hpp"
#include "apt/vec.hpp"
#include "apt/virtual_vec.hpp"
#include "behaviors.hpp"
#include "catch2/catch.hpp"

using namespace apt;

SCENARIO("vec-vec algebra", "[apt]") {
  Vec<int,3> v1 ( 4, 9, 12 );
  Vec<int,3> v2 ( 1, 3, 2 );
  WHEN("v1 + v2") {
    REQUIRE_FALSE(is<bhv::lvec>(v1 + v2));

    THEN("assignalbe to Vec") {
      Vec<int,3> v_sum = v1 + v2;
      REQUIRE( std::get<0>(v_sum) == 5 );
      REQUIRE( std::get<1>(v_sum) == 12 );
      REQUIRE( std::get<2>(v_sum) == 14 );
    }
  }

  WHEN("v1 - v2") {
    REQUIRE_FALSE(is<bhv::lvec>(v1 - v2));

    THEN("assignalbe to Vec") {
      Vec<int,3> v_sub = v1 - v2;
      REQUIRE( std::get<0>(v_sub) == 3 );
      REQUIRE( std::get<1>(v_sub) == 6 );
      REQUIRE( std::get<2>(v_sub) == 10 );
    }
  }

  WHEN("v1 * v2") {
    REQUIRE_FALSE(is<bhv::lvec>(v1 * v2));

    THEN("assignalbe to Vec") {
      Vec<int,3> v_mul = v1 * v2;
      REQUIRE( std::get<0>(v_mul) == 4 );
      REQUIRE( std::get<1>(v_mul) == 27 );
      REQUIRE( std::get<2>(v_mul) == 24 );
    }
  }

  WHEN("v1 / v2") {
    REQUIRE_FALSE(is<bhv::lvec>(v1 / v2));

    THEN("assignalbe to Vec") {
      Vec<int,3> v_div = v1 / v2;
      REQUIRE( std::get<0>(v_div) == 4 );
      REQUIRE( std::get<1>(v_div) == 3 );
      REQUIRE( std::get<2>(v_div) == 6 );
    }
  }

}

SCENARIO("vec-sca algebra", "[apt]") {
  Vec<int,3> v ( 3, 9, 12 );
  int s = 3;
  WHEN("v + s") {
    REQUIRE_FALSE(is<bhv::lvec>(v + s));

    THEN("assignalbe to Vec") {
      Vec<int,3> v_sum = v + s;
      REQUIRE( std::get<0>(v_sum) == 6 );
      REQUIRE( std::get<1>(v_sum) == 12 );
      REQUIRE( std::get<2>(v_sum) == 15 );
    }
  }

  WHEN("v - s") {
    REQUIRE_FALSE(is<bhv::lvec>(v - s));

    THEN("assignalbe to Vec") {
      Vec<int,3> v_sub = v - s;
      REQUIRE( std::get<0>(v_sub) == 0 );
      REQUIRE( std::get<1>(v_sub) == 6 );
      REQUIRE( std::get<2>(v_sub) == 9 );
    }
  }

  WHEN("v * s") {
    REQUIRE_FALSE(is<bhv::lvec>(v * s));

    THEN("assignalbe to Vec") {
      Vec<int,3> v_mul = v * s;
      REQUIRE( std::get<0>(v_mul) == 9 );
      REQUIRE( std::get<1>(v_mul) == 27 );
      REQUIRE( std::get<2>(v_mul) == 36 );
    }
  }

  WHEN("v / s") {
    REQUIRE_FALSE(is<bhv::lvec>(v / s));

    THEN("assignalbe to Vec") {
      Vec<int,3> v_div = v / s;
      REQUIRE( std::get<0>(v_div) == 1 );
      REQUIRE( std::get<1>(v_div) == 3 );
      REQUIRE( std::get<2>(v_div) == 4 );
    }
  }

}

SCENARIO("vec algebra with assign", "[apt]") {
  Vec<int,3> v1 ( 7, 70, 700 );
  Vec<int,3> v2 ( 1, 3, 2 );
  int s = 7;
  WHEN("v1 += v2") {
    // TODO fix these
    // REQUIRE(is<bhv::lvec>(v1 += v2));

    v1 += v2;
    REQUIRE( std::get<0>(v1) == 8 );
    REQUIRE( std::get<1>(v1) == 73 );
    REQUIRE( std::get<2>(v1) == 702 );

    std::get<0>(v1 += v2) = 17;
    REQUIRE( std::get<0>(v1) == 17 );
    REQUIRE( std::get<1>(v1) == 76 );
    REQUIRE( std::get<2>(v1) == 704 );
  }

  WHEN("v1 -= v2") {
    // REQUIRE(is<bhv::lvec>(v1 -= v2));

    v1 -= v2;
    REQUIRE( std::get<0>(v1) == 6 );
    REQUIRE( std::get<1>(v1) == 67 );
    REQUIRE( std::get<2>(v1) == 698 );

    std::get<0>(v1 -= v2) = 17;
    REQUIRE( std::get<0>(v1) == 17 );
    REQUIRE( std::get<1>(v1) == 64 );
    REQUIRE( std::get<2>(v1) == 696 );
  }

  WHEN("v1 *= s") {
    // REQUIRE(is<bhv::lvec>(v1 *= s));

    v1 *= s;
    REQUIRE( std::get<0>(v1) == 49 );
    REQUIRE( std::get<1>(v1) == 490 );
    REQUIRE( std::get<2>(v1) == 4900 );
  }

  WHEN("v1 /= s") {
    // REQUIRE(is<bhv::lvec>(v1 /= s));

    v1 /= s;
    REQUIRE( std::get<0>(v1) == 1 );
    REQUIRE( std::get<1>(v1) == 10 );
    REQUIRE( std::get<2>(v1) == 100 );
  }
}

SCENARIO("vec cross product", "[apt]") {
  Vec<int,3> v1 ( 7, 9, -5 );
  Vec<int,3> v2 ( -3, 2, -13 );
  auto v_crs = cross(v1, v2);
  REQUIRE_FALSE( is<bhv::lvec>(v_crs) );
  REQUIRE( std::get<0>(v_crs) == -107 );
  REQUIRE( std::get<1>(v_crs) == 106 );
  REQUIRE( std::get<2>(v_crs) == 41 );
}

SCENARIO("vec inner product", "[apt]") {
  Vec<int,3> v1 ( 7, 9, -5 );
  Vec<int,3> v2 ( -3, 2, -13 );

  REQUIRE( dot(v1,v2) == 62 );
  REQUIRE( sqabs(v1) == 155 );
  REQUIRE( abs(v2) == Approx(std::sqrt(182)) );
}

