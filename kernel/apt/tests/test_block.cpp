#include "testfw/testfw.hpp"
#include "apt/block.hpp"

using namespace apt;

SCENARIO("Block", "[apt]") {
  Index<3> beg { 2, 3, 4};
  Index<3> end { 13, 15, 17};
  Block block(beg, end);
  auto itr = block.begin();
  REQUIRE( block.end() == end[2] );

  for ( int k = beg[2]; k < end[2]; ++k )
    for ( int j = beg[1]; j < end[1]; ++j )
      for ( int i = beg[0]; i < end[0]; ++i ) {
        CAPTURE(*itr, i, j, k);
        REQUIRE( *(itr++) == Index<3>{ i, j, k } );
      }
}

TEMPLATE_TEST_CASE("Test empty blocks", "[apt]"
                   , (aio::IndexType<0,0,0>)
                   , (aio::IndexType<1,0,0>)
                   , (aio::IndexType<0,1,0>)
                   , (aio::IndexType<0,0,1>)
                   , (aio::IndexType<0,1,1>)
                   , (aio::IndexType<1,0,1>)
                   , (aio::IndexType<1,1,0>)
                   ) {
  constexpr int DGrid = 3;
  Block<DGrid> block( {}, TestType::get() );
  // block.begin() should equal to block.end()
  REQUIRE_FALSE( ( block.begin() != block.end() ) );
}
