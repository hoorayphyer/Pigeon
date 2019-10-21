#include "testfw/testfw.hpp"
#include "field/mesh.hpp"

using namespace field;

SCENARIO("Test LinearizedIndex", "[field][mesh]") {
  constexpr int D = 2;
  apt::Range r0 ( 0, 128, {2,3} );
  apt::Range r1 ( 0, 256, {1,2} );
  Mesh<D> mesh( {r0,r1} );
  REQUIRE(mesh.stride()[0] == 1);
  REQUIRE(mesh.stride()[1] == 128 + 2 + 3);
  REQUIRE(mesh.stride()[2] == (128 + 2 + 3) * (256 + 1 + 2));

  aio::unif_int<unsigned int> dist(0, 255);

  int N = 100*1000;
  while( N-- ) {
    apt::Index<D> data { dist(), dist() }; // NOTE the index can go out of bounds, but it is fine here.
    auto index = mesh.linear_index(data);
    REQUIRE( index == (data[0] + r0.margin()[LFT]) * mesh.stride()[0] +  (data[1] + r1.margin()[LFT]) * mesh.stride()[1] );
  }

}
