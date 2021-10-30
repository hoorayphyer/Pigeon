#include <limits>

#include "particle/state.hpp"
#include "testfw/testfw.hpp"

using namespace particle;

using T = long long;

struct State : public StateExpression<State, T> {
 private:
  T _state = 0;

 public:
  constexpr T& state() noexcept { return _state; }
  constexpr T state() const noexcept { return _state; }
};

SCENARIO("setting species", "[particle]") {
  State s;
  REQUIRE(s.get<species>() == species::unknown);

  s.set(species::electron);
  REQUIRE(s.get<species>() == species::electron);
  s.set(species::positron);
  REQUIRE(s.get<species>() == species::positron);
  s.set(species::ion);
  REQUIRE(s.get<species>() == species::ion);
  s.set(species::photon);
  REQUIRE(s.get<species>() == species::photon);
}

SCENARIO("setting individual flags", "[particle]") {
  {
    State s;
    REQUIRE_FALSE(s.is(flag::exist));
  }

  auto f = [](flag X) {
    State s;
    REQUIRE_FALSE(s.is(X));
    s.set(X);
    REQUIRE(s.is(X));
    s.reset(X);
    REQUIRE_FALSE(s.is(X));
  };

  static_assert(2 * sizeof(flagbits) == 16);
  for (std::underlying_type_t<flag> i = 0; i < 2 * sizeof(flagbits); ++i) {
    f(static_cast<flag>(i));
  }
}

TEMPLATE_TEST_CASE(
    "Test assigning flags", "[particle]", (std::integral_constant<int, 0>),
    (std::integral_constant<int, 1>), (std::integral_constant<int, 2>),
    (std::integral_constant<int, 3>), (std::integral_constant<int, 4>),
    (std::integral_constant<int, 5>), (std::integral_constant<int, 6>),
    (std::integral_constant<int, 7>), (std::integral_constant<int, 8>),
    (std::integral_constant<int, 9>), (std::integral_constant<int, 10>),
    (std::integral_constant<int, 11>), (std::integral_constant<int, 12>),
    (std::integral_constant<int, 13>), (std::integral_constant<int, 14>),
    (std::integral_constant<int, 15>)) {
  constexpr auto FLAG = static_cast<flag>(TestType::value);

  State s;
  CHECK_FALSE(s.is(FLAG));
  s.assign<FLAG>(true);
  CHECK(s.is(FLAG));
  s.assign<FLAG>(false);
  CHECK_FALSE(s.is(FLAG));
}

SCENARIO("Test setting individual bits of flagbits", "[particle]") {
  static_assert(2 * sizeof(flagbits) == 16);
  for (std::underlying_type_t<flag> i = 0; i < 2 * sizeof(flagbits); ++i) {
    flagbits x{};
    x[static_cast<flag>(i)] = true;
    CHECK(x.to_ulong() == (1u << i));
    x[static_cast<flag>(i)] = false;
    CHECK(x.to_ulong() == 0);
  }
}

TEMPLATE_TEST_CASE("Test setting migrcode", "[particle]",
                   (std::integral_constant<int, 1>),
                   (std::integral_constant<int, 2>),
                   (std::integral_constant<int, 3>)) {
  constexpr int D = TestType::value;
  const int max = apt::pow3(D);

  State s;
  for (int i = 0; i < max; ++i) {
    s.set<migrcode, D>(i);
    CHECK(s.get<migrcode, D>() == i);
    s.reset<migrcode>();
    CHECK(s.get<migrcode, D>() == (apt::pow3(D) - 1) / 2);
  }
}

TEMPLATE_TEST_CASE("Test set and get attributes", "[particle]", flagbits,
                   pid  // NOTE this may take a few seconds
) {
  State s;
  std::size_t max = (1u << layout::size<TestType>());

  for (std::size_t i = 0; i < max; ++i) {
    TestType attr(i);
    s.set(attr);
    CAPTURE(i);
    CHECK(s.get<TestType>() == i);
    s.reset<TestType>();
    CHECK(s.get<TestType>() == 0);
  }
}

SCENARIO("Test setting mixed attributions", "[particle]") {
  State s;
  s.set(species::electron, flagbits(147), flag::secondary, pid(7468));
  s.set<migrcode, 2>(17);
  CHECK(s.get<species>() == species::electron);
  CHECK(s.get<flagbits>() == (147 | (1 << static_cast<int>(flag::secondary))));
  CHECK(s.get<pid>() == 7468);
  CHECK(s.get<migrcode, 2>() == 17);
}
