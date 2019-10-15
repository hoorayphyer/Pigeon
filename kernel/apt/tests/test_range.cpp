#include "testfw/testfw.hpp"
#include "apt/range.hpp"

using namespace apt;

SCENARIO("Test single range", "[apt]") {
  int begin = 13, end = 56;
  apt::pair<int> margin = {3, 2};
  Range r ( begin, end, margin );
  REQUIRE( r.begin() == begin );
  REQUIRE( r.end() == end );
  REQUIRE( r.margin()[LFT] == margin[LFT] );
  REQUIRE( r.margin()[RGT] == margin[RGT] );
  REQUIRE( r.size() == end - begin );
  REQUIRE( r.far_begin() == begin - margin[LFT] );
  REQUIRE( r.far_end() == end + margin[RGT] );
  REQUIRE( r.full_size() == end + margin[RGT] - begin + margin[LFT] );
}

SCENARIO("Test range array", "[apt]") {
  constexpr int D = 2;
  Range r1 ( 13, 78, {2,3} );
  Range r2 ( 56, 397, {4,1} );
  {
    auto x = range::begin<D>({r1,r2});
    REQUIRE(x[0] == r1.begin());
    REQUIRE(x[1] == r2.begin());
  }
  {
    auto x = range::end<D>({r1,r2});
    REQUIRE(x[0] == r1.end());
    REQUIRE(x[1] == r2.end());
  }
  {
    auto x = range::size<D>({r1,r2});
    REQUIRE(x[0] == r1.size());
    REQUIRE(x[1] == r2.size());
  }
  {
    auto x = range::far_begin<D>({r1,r2});
    REQUIRE(x[0] == r1.far_begin());
    REQUIRE(x[1] == r2.far_begin());
  }
  {
    auto x = range::far_end<D>({r1,r2});
    REQUIRE(x[0] == r1.far_end());
    REQUIRE(x[1] == r2.far_end());
  }
  {
    auto x = range::full_size<D>({r1,r2});
    REQUIRE(x[0] == r1.full_size());
    REQUIRE(x[1] == r2.full_size());
  }
  REQUIRE_FALSE( range::is_margin_uniform<D>( {r1,r2} ) );
  r1.margin() = {7,7};
  r2.margin() = {7,7};
  REQUIRE( range::is_margin_uniform<D>( {r1,r2} ) );

}
