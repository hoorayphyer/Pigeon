#include "testfw/testfw.hpp"
#include "timer/timer.hpp"
#include <unistd.h>

using namespace tmr;

SCENARIO("Test Timer", "[timer]") {
  SECTION("Test lapse") {
    Timestamp stamp1;
    sleep(1);
    REQUIRE(stamp1.lapse() == Approx(1000).epsilon(0.01));
  }
  SECTION("Test operator-") {
    Timestamp stamp1;
    sleep(1);
    Timestamp stamp2;
    REQUIRE(stamp2 - stamp1 == Approx(1000).epsilon(0.01));
  }
}
