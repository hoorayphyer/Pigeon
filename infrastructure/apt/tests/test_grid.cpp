#include "apt/grid.hpp"
#include "testfw/testfw.hpp"

using namespace apt;
using T = float;

SCENARIO("grid csba", "[apt]") {
  Grid1D<T> g(0.3, 1.2, 17);
  CHECK(g.csba(0.3) == 0);
  CHECK(g.csba(0.3 + 0.91 / 17) == 1);
  CHECK(g.csba(0.3 - 0.89 / 17) == -1);

  for (int i = -1024; i < 1025; ++i) {
    T x = g.absc(i, 0.5);
    int j = g.csba(x);
    CAPTURE(i, j);
    CHECK(i == j);
  }
}
