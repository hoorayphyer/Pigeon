#include "testfw/testfw.hpp"
#include "apt/block.hpp"

using namespace apt;

SCENARIO("Block", "[apt]") {
  int N[3] = { 13, 15, 17 };
  Block block{ Index<3>{ N[0], N[1], N[2] }};
  auto itr = block.begin();
  auto end = block.end();
  REQUIRE_FALSE( end != Index<3>{ 0, 0, N[2] } );

  for ( int k = 0; k < N[2]; ++k )
    for ( int j = 0; j < N[1]; ++j )
      for ( int i = 0; i < N[0]; ++i ) {
        CAPTURE(*itr, i, j, k);
        REQUIRE( *(itr++) == Index<3>{ i, j, k } );
      }
}

TEMPLATE_TEST_CASE("Testing deposition in 2D against alternative implementation BrutalForce_dJ_Field", "[field][mpi]"
                   , (aio::IndexType<0,0,0>)
                   , (aio::IndexType<1,0,0>)
                   , (aio::IndexType<0,1,0>)
                   , (aio::IndexType<0,0,1>)
                   , (aio::IndexType<0,1,1>)
                   , (aio::IndexType<1,0,1>)
                   , (aio::IndexType<1,1,0>)
                   ) {
  constexpr int DGrid = 3;
  Block<DGrid> block( TestType::get() );
  REQUIRE_FALSE( block.begin() != block.end());
}
