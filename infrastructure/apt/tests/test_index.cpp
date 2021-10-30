#include "apt/index.hpp"
#include "testfw/testfw.hpp"

using namespace apt;

SCENARIO("Test Longidx constructors", "[apt]") {
  Longidx x;
  REQUIRE(x == 0);
  REQUIRE(x.dir() == 0);

  Longidx y(1, -6);
  REQUIRE(y == -6);
  REQUIRE(y.dir() == 1);
}

SCENARIO("Test Longidx arithmetics", "[apt]") {
  SECTION("unary") {
    Longidx x(5, 6);
    REQUIRE(++x == 7);
    REQUIRE(x == 7);
    REQUIRE(--x == 6);
    REQUIRE(x == 6);
    REQUIRE(x++ == 6);
    REQUIRE(x == 7);
    REQUIRE(x-- == 7);
    REQUIRE(x == 6);
  }

  SECTION("comparison") {
    REQUIRE(Longidx(5, 2) < 3);
    REQUIRE(Longidx(5, 2) <= 2);
    REQUIRE(Longidx(5, 2) > 1);
    REQUIRE(Longidx(5, 2) >= 2);
    REQUIRE(Longidx(5, 2) != 9);
    REQUIRE(Longidx(5, 2) == 2);
  }

  SECTION("binary, with Index<D>") {
    constexpr int D = 3;
    Index<D> b{2, 13, -91};
    Longidx a(1, 6);
    REQUIRE(a + b == Index<D>{2, 19, -91});
    REQUIRE(b + a == Index<D>{2, 19, -91});
    REQUIRE(b - a == Index<D>{2, 7, -91});
    REQUIRE(a - b == Index<D>{-2, -7, 91});
    Longidx c(2, -8);
    REQUIRE(c + b == Index<D>{2, 13, -99});
    REQUIRE(b + c == Index<D>{2, 13, -99});
    REQUIRE(b - c == Index<D>{2, 13, -83});
    REQUIRE(c - b == Index<D>{-2, -13, 83});
  }

  SECTION("binary, with int") {
    Longidx a(1, 147);
    REQUIRE(a + 5 == 152);
    REQUIRE(a - 29 == 118);
    REQUIRE(a * 14 == 147 * 14);
    REQUIRE(a / 37 == 147 / 37);

    REQUIRE(5 + a == 152);
    REQUIRE(29 - a == -118);
    REQUIRE(14 * a == 147 * 14);
  }
}

SCENARIO("Test iteration with Longidx", "[apt]") {
  int i = 0;
  for (Longidx n{1, 6}; n < 18; ++n) {
    REQUIRE(n == 6 + i);
    ++i;
  }
}
