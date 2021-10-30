#include <unistd.h>

#include "testfw/testfw.hpp"
#include "timer/timer.hpp"

using namespace tmr;

SCENARIO("Test Timer", "[timer]") {
  SECTION("Test lapse") {
    Timestamp stamp1;
    sleep(1);
    REQUIRE(stamp1.lapse().in_units_of("s").val() == Approx(1).epsilon(0.01));
    REQUIRE(stamp1.lapse().in_units_of("ms").val() ==
            Approx(1000).epsilon(0.01));
    REQUIRE(stamp1.lapse().in_units_of("us").val() ==
            Approx(1000 * 1000).epsilon(0.01));
    REQUIRE(stamp1.lapse().in_units_of("ns").val() ==
            Approx(1000 * 1000 * 1000).epsilon(0.01));
  }
  SECTION("Test operator-") {
    Timestamp stamp1;
    sleep(1);
    Timestamp stamp2;
    REQUIRE((stamp2 - stamp1).in_units_of("ms").val() ==
            Approx(1000).epsilon(0.01));
  }
}
