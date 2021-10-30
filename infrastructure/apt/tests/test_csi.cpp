#include "apt/csi.hpp"
#include "testfw/testfw.hpp"

using namespace apt;

// comma-separated-integer
TEMPLATE_TEST_CASE("Test csi", "[apt]", int, long, long long) {
  using T = TestType;
  T x{};
  SECTION("beginning zeros") {
    x = 1024;
    CHECK(csi(x) == "1,024");
    CHECK_FALSE(csi(x) == "1,24");

    x = 1004;
    CHECK(csi(x) == "1,004");
    CHECK_FALSE(csi(x) == "1,4");

    x = 1000;
    CHECK(csi(x) == "1,000");
    CHECK_FALSE(csi(x) == "1,0");
    CHECK_FALSE(csi(x) == "1,");
  }

  SECTION("positive number") {
    x = 7;
    CHECK(csi(x) == "7");
    CHECK_FALSE(csi(x) == ",7");
    x = 47;
    CHECK(csi(x) == "47");
    x = 147;
    CHECK(csi(x) == "147");

    x = 6147;
    CHECK(csi(x) == "6,147");
    x = 76147;
    CHECK(csi(x) == "76,147");
    x = 376147;
    CHECK(csi(x) == "376,147");

    x = 9376147;
    CHECK(csi(x) == "9,376,147");
    x = 59376147;
    CHECK(csi(x) == "59,376,147");
    x = 859376147;
    CHECK(csi(x) == "859,376,147");

    if (!std::is_same_v<T, int>) {
      x = 2859376147;
      CHECK(csi(x) == "2,859,376,147");
      x = 22859376147;
      CHECK(csi(x) == "22,859,376,147");
      x = 122859376147;
      CHECK(csi(x) == "122,859,376,147");

      x = 4122859376147;
      CHECK(csi(x) == "4,122,859,376,147");
      x = 94122859376147;
      CHECK(csi(x) == "94,122,859,376,147");
      x = 794122859376147;
      CHECK(csi(x) == "794,122,859,376,147");

      x = 5794122859376147;
      CHECK(csi(x) == "5,794,122,859,376,147");
    }
  }

  SECTION("negative number") {
    x = -0;
    CHECK(csi(x) == "0");
    CHECK_FALSE(csi(x) == "-0");
    x = -7;
    CHECK(csi(x) == "-7");
    CHECK_FALSE(csi(x) == "-,7");
    CHECK_FALSE(csi(x) == ",-7");
    x = -47;
    CHECK(csi(x) == "-47");
    x = -147;
    CHECK(csi(x) == "-147");

    x = -6147;
    CHECK(csi(x) == "-6,147");
    x = -76147;
    CHECK(csi(x) == "-76,147");
    x = -376147;
    CHECK(csi(x) == "-376,147");

    x = -9376147;
    CHECK(csi(x) == "-9,376,147");
    x = -59376147;
    CHECK(csi(x) == "-59,376,147");
    x = -859376147;
    CHECK(csi(x) == "-859,376,147");

    if (!std::is_same_v<T, int>) {
      x = -2859376147;
      CHECK(csi(x) == "-2,859,376,147");
      x = -22859376147;
      CHECK(csi(x) == "-22,859,376,147");
      x = -122859376147;
      CHECK(csi(x) == "-122,859,376,147");

      x = -4122859376147;
      CHECK(csi(x) == "-4,122,859,376,147");
      x = -94122859376147;
      CHECK(csi(x) == "-94,122,859,376,147");
      x = -794122859376147;
      CHECK(csi(x) == "-794,122,859,376,147");

      x = -5794122859376147;
      CHECK(csi(x) == "-5,794,122,859,376,147");
    }
  }
}
