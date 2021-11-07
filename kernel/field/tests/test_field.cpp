#include "field/field.hpp"
#include "testfw/testfw.hpp"

using namespace field;

SCENARIO("Test setting offset") {
  Mesh<2> mesh({{{0, 32, 1}, {0, 32, 1}}});
  Field<double, 3, 2> f(mesh);
  apt::array<offset_t, 2> offset{MIDWAY, MIDWAY};
  for (int i = 0; i < 3; ++i) {
    f.set_offset(i, offset);
    REQUIRE(f[i].offset() == offset);
  }

  REQUIRE((!MIDWAY) == INSITU);
  REQUIRE((!INSITU) == MIDWAY);
}
