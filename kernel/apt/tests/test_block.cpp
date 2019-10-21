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

SCENARIO("Test empty blocks", "[apt]") {
  constexpr int D = 3;
  auto is_empty =
    []( const Block<D>& b ) noexcept {
      return !( b.begin() != b.end() );
    };
  REQUIRE_FALSE( is_empty({ {0,0,0}, {1,1,1} } ) );
  REQUIRE( is_empty({ {1,0,0}, {1,1,1} } ) );
  REQUIRE( is_empty({ {2,0,0}, {1,1,1} } ) );
  REQUIRE( is_empty({ {0,1,0}, {1,1,1} } ) );
  REQUIRE( is_empty({ {0,2,0}, {1,1,1} } ) );
  REQUIRE( is_empty({ {0,0,1}, {1,1,1} } ) );
  REQUIRE( is_empty({ {0,0,2}, {1,1,1} } ) );

  WHEN("loop over empty block") {
    Block<D> block( { 12, 23, 34}, { 12, 26, 37 });
    int n = 0;
    for ( const auto& I : block ) ++n;
    REQUIRE( 0 == n );
  }
}
