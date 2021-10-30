#include "apt/array.hpp"
#include "apt/foreach.hpp"
#include "testfw/testfw.hpp"

using namespace apt;

SCENARIO("lambda f takes auto&", "[apt]") {
  auto f = [](auto& a) { a *= 2; };
  array<int, 3> x{147, 147, 147};

  foreach
    <0, 3>(f, x);
  for (int i = 0; i < 3; ++i) {
    REQUIRE(x[i] == 147 * 2);
  }
}

SCENARIO("lambda f takes auto", "[apt]") {
  auto f = [](auto a) { a *= 2; };
  array<int, 3> x{147, 147, 147};

  foreach
    <0, 3>(f, x);
  for (int i = 0; i < 3; ++i) {
    REQUIRE(x[i] == 147);
  }
}

SCENARIO("lambda f takes auto&&", "[apt]") {
  WHEN("the argument is not modified") {
    auto f = [](auto&& a) { a *= 2; };
    array<int, 3> x{147, 147, 147};

    foreach
      <0, 3>(f, x);
    for (int i = 0; i < 3; ++i) {
      REQUIRE(x[i] == 147 * 2);
    }
  }

  WHEN("the argument is modified") {
    auto f = [](auto&& a) { a * 2; };
    array<int, 3> x{147, 147, 147};

    foreach
      <0, 3>(f, x);
    for (int i = 0; i < 3; ++i) {
      REQUIRE(x[i] == 147);
    }
  }
}
