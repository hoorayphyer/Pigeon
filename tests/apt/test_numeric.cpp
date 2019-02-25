#include "apt/numeric.hpp"
#include "apt/vec.hpp"
#include "apt/virtual_vec.hpp"
#include "behaviors.hpp"
#include "catch2/catch.hpp"

using namespace apt;

SCENARIO("vec-vec algebra", "[apt][vec]") {
  Vec<int,3> v1 ( 4, 9, 12 );
  Vec<int,3> v2 ( 1, 3, 2 );
  WHEN("v1 + v2") {
    REQUIRE_FALSE(is<bhv::lvec>(v1 + v2));

    THEN("assignalbe to Vec") {
      Vec<int,3> v_sum = v1 + v2;
      REQUIRE( v_sum[0] == 5 );
      REQUIRE( v_sum[1] == 12 );
      REQUIRE( v_sum[2] == 14 );
    }
  }

  WHEN("v1 - v2") {
    REQUIRE_FALSE(is<bhv::lvec>(v1 - v2));

    THEN("assignalbe to Vec") {
      Vec<int,3> v_sub = v1 - v2;
      REQUIRE( v_sub[0] == 3 );
      REQUIRE( v_sub[1] == 6 );
      REQUIRE( v_sub[2] == 10 );
    }
  }

  WHEN("v1 * v2") {
    REQUIRE_FALSE(is<bhv::lvec>(v1 * v2));

    THEN("assignalbe to Vec") {
      Vec<int,3> v_mul = v1 * v2;
      REQUIRE( v_mul[0] == 4 );
      REQUIRE( v_mul[1] == 27 );
      REQUIRE( v_mul[2] == 24 );
    }
  }

  WHEN("v1 / v2") {
    REQUIRE_FALSE(is<bhv::lvec>(v1 / v2));

    THEN("assignalbe to Vec") {
      Vec<int,3> v_div = v1 / v2;
      REQUIRE( v_div[0] == 4 );
      REQUIRE( v_div[1] == 3 );
      REQUIRE( v_div[2] == 6 );
    }
  }

}

SCENARIO("vec-sca algebra", "[apt][vec]") {
  Vec<int,3> v ( 3, 9, 12 );
  int s = 3;
  WHEN("v + s") {
    REQUIRE_FALSE(is<bhv::lvec>(v + s));

    THEN("assignalbe to Vec") {
      Vec<int,3> v_sum = v + s;
      REQUIRE( v_sum[0] == 6 );
      REQUIRE( v_sum[1] == 12 );
      REQUIRE( v_sum[2] == 15 );
    }
  }

  WHEN("v - s") {
    REQUIRE_FALSE(is<bhv::lvec>(v - s));

    THEN("assignalbe to Vec") {
      Vec<int,3> v_sub = v - s;
      REQUIRE( v_sub[0] == 0 );
      REQUIRE( v_sub[1] == 6 );
      REQUIRE( v_sub[2] == 9 );
    }
  }

  WHEN("v * s") {
    REQUIRE_FALSE(is<bhv::lvec>(v * s));

    THEN("assignalbe to Vec") {
      Vec<int,3> v_mul = v * s;
      REQUIRE( v_mul[0] == 9 );
      REQUIRE( v_mul[1] == 27 );
      REQUIRE( v_mul[2] == 36 );
    }
  }

  WHEN("v / s") {
    REQUIRE_FALSE(is<bhv::lvec>(v / s));

    THEN("assignalbe to Vec") {
      Vec<int,3> v_div = v / s;
      REQUIRE( v_div[0] == 1 );
      REQUIRE( v_div[1] == 3 );
      REQUIRE( v_div[2] == 4 );
    }
  }

}

SCENARIO("vec algebra with assign", "[apt][vec]") {
  Vec<int,3> v1 ( 7, 70, 700 );
  Vec<int,3> v2 ( 1, 3, 2 );
  int s = 7;
  WHEN("v1 += v2") {
    // TODO fix these
    // REQUIRE(is<bhv::lvec>(v1 += v2));

    v1 += v2;
    REQUIRE( v1[0] == 8 );
    REQUIRE( v1[1] == 73 );
    REQUIRE( v1[2] == 702 );

    std::get<0>(v1 += v2) = 17;
    REQUIRE( v1[0] == 17 );
    REQUIRE( v1[1] == 76 );
    REQUIRE( v1[2] == 704 );
  }

  WHEN("v1 -= v2") {
    // REQUIRE(is<bhv::lvec>(v1 -= v2));

    v1 -= v2;
    REQUIRE( v1[0] == 6 );
    REQUIRE( v1[1] == 67 );
    REQUIRE( v1[2] == 698 );

    std::get<0>(v1 -= v2) = 17;
    REQUIRE( v1[0] == 17 );
    REQUIRE( v1[1] == 64 );
    REQUIRE( v1[2] == 696 );
  }

  WHEN("v1 *= s") {
    // REQUIRE(is<bhv::lvec>(v1 *= s));

    v1 *= s;
    REQUIRE( v1[0] == 49 );
    REQUIRE( v1[1] == 490 );
    REQUIRE( v1[2] == 4900 );
  }

  WHEN("v1 /= s") {
    // REQUIRE(is<bhv::lvec>(v1 /= s));

    v1 /= s;
    REQUIRE( v1[0] == 1 );
    REQUIRE( v1[1] == 10 );
    REQUIRE( v1[2] == 100 );
  }
}

SCENARIO("vec cross product", "[apt][vec]") {
  Vec<int,3> v1 ( 7, 9, -5 );
  Vec<int,3> v2 ( -3, 2, -13 );
  auto v_crs = cross(v1, v2);
  REQUIRE_FALSE( is<bhv::lvec>(v_crs) );
  REQUIRE( v_crs[0] == -107 );
  REQUIRE( v_crs[1] == 106 );
  REQUIRE( v_crs[2] == 41 );
}

SCENARIO("vec inner product", "[apt][vec]") {
  Vec<int,3> v1 ( 7, 9, -5 );
  Vec<int,3> v2 ( -3, 2, -13 );

  REQUIRE( dot(v1,v2) == 62 );
  REQUIRE( sqabs(v1) == 155 );
  REQUIRE( abs(v2) == Approx(std::sqrt(182)) );
}

