#include "apt/range.hpp"
#include "testfw/testfw.hpp"

using namespace apt;

SCENARIO("Test Range", "[apt]") {
  int begin = 13, end = 56;
  apt::pair<int> margin = {3, 2};
  Range r(begin, end, margin);
  REQUIRE(r.begin() == begin);
  REQUIRE(r.end() == end);
  REQUIRE(r.margin()[LFT] == margin[LFT]);
  REQUIRE(r.margin()[RGT] == margin[RGT]);
  REQUIRE(r.size() == end - begin);
  REQUIRE(r.far_begin() == begin - margin[LFT]);
  REQUIRE(r.far_end() == end + margin[RGT]);
  REQUIRE(r.full_size() == end + margin[RGT] - begin + margin[LFT]);
}

SCENARIO("Test Range comparison", "[apt]") {
  int begin = 13, end = 56;
  apt::pair<int> margin = {3, 2};
  Range r1(begin, end, margin);
  Range r2 = {r1};
  REQUIRE((r1 == r2));
  REQUIRE_FALSE((r1 != r2));
  r2.begin() = 14;
  REQUIRE((r1 != r2));
  r2.begin() = r1.begin();
  REQUIRE((r1 == r2));
  r2.end() = 55;
  REQUIRE((r1 != r2));
  r2.end() = 56;
  REQUIRE((r1 == r2));
  r2.margin()[LFT] = 2;
  REQUIRE((r1 != r2));
  r2.margin()[LFT] = 3;
  REQUIRE((r1 == r2));
  r2.margin()[RGT] = 3;
  REQUIRE((r1 != r2));
  r2.margin()[RGT] = 2;
  REQUIRE((r1 == r2));
}

SCENARIO("Test array of range ", "[apt]") {
  constexpr int D = 2;
  Range r1(13, 78, {2, 3});
  Range r2(56, 397, {4, 1});
  array<Range, D> rgs{r1, r2};
  {
    auto x = range::begin<D>(rgs);
    REQUIRE(x[0] == r1.begin());
    REQUIRE(x[1] == r2.begin());

    REQUIRE(range::begin(rgs, 0) == r1.begin());
    REQUIRE(range::begin(rgs, 1) == r2.begin());
  }
  {
    auto x = range::end<D>(rgs);
    REQUIRE(x[0] == r1.end());
    REQUIRE(x[1] == r2.end());

    REQUIRE(range::end(rgs, 0) == r1.end());
    REQUIRE(range::end(rgs, 1) == r2.end());
  }
  {
    auto x = range::margin<D>(rgs);
    REQUIRE(x[0][LFT] == r1.margin()[LFT]);
    REQUIRE(x[0][RGT] == r1.margin()[RGT]);
    REQUIRE(x[1][LFT] == r2.margin()[LFT]);
    REQUIRE(x[1][RGT] == r2.margin()[RGT]);

    REQUIRE(range::margin(rgs, 0)[LFT] == r1.margin()[LFT]);
    REQUIRE(range::margin(rgs, 0)[RGT] == r1.margin()[RGT]);
    REQUIRE(range::margin(rgs, 1)[LFT] == r2.margin()[LFT]);
    REQUIRE(range::margin(rgs, 1)[RGT] == r2.margin()[RGT]);
  }
  {
    auto x = range::size<D>(rgs);
    REQUIRE(x[0] == r1.size());
    REQUIRE(x[1] == r2.size());

    REQUIRE(range::size(rgs, 0) == r1.size());
    REQUIRE(range::size(rgs, 1) == r2.size());
  }
  {
    auto x = range::far_begin<D>(rgs);
    REQUIRE(x[0] == r1.far_begin());
    REQUIRE(x[1] == r2.far_begin());

    REQUIRE(range::far_begin(rgs, 0) == r1.far_begin());
    REQUIRE(range::far_begin(rgs, 1) == r2.far_begin());
  }
  {
    auto x = range::far_end<D>(rgs);
    REQUIRE(x[0] == r1.far_end());
    REQUIRE(x[1] == r2.far_end());

    REQUIRE(range::far_end(rgs, 0) == r1.far_end());
    REQUIRE(range::far_end(rgs, 1) == r2.far_end());
  }
  {
    auto x = range::full_size<D>(rgs);
    REQUIRE(x[0] == r1.full_size());
    REQUIRE(x[1] == r2.full_size());

    REQUIRE(range::full_size(rgs, 0) == r1.full_size());
    REQUIRE(range::full_size(rgs, 1) == r2.full_size());
  }
  REQUIRE_FALSE(range::is_margin_uniform<D>(rgs));
  rgs[0].margin() = {7, 7};
  rgs[1].margin() = {7, 7};
  REQUIRE(range::is_margin_uniform<D>(rgs));
}
