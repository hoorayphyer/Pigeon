#include "apt/pair.hpp"
#include "testfw/testfw.hpp"

using namespace apt;

SCENARIO("test pair", "[apt]") {
  pair<int> x{6, 8};
  REQUIRE(x[0] == 6);
  REQUIRE(x[LFT] == 6);
  REQUIRE(x[1] == 8);
  REQUIRE(x[RGT] == 8);
  x[0]++;
  REQUIRE(x[0] == 7);
}

SCENARIO("test tie", "[apt]") {
  pair<int> x{5, 3};
  int a = 0, b = 0;
  tie(a, b) = x;
  REQUIRE(a == 5);
  REQUIRE(b == 3);
}
