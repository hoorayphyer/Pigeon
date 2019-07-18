#include "testfw/testfw.hpp"
#include "particle/state.hpp"
#include <limits>

using namespace particle;

using T = long long;

struct State : public StateExpression<State,T> {
private:
  T _state = 0;
public:
  constexpr T& state() noexcept { return _state; }
  constexpr T state() const noexcept { return _state; }
};

TEMPLATE_TEST_CASE("Test setbits and getbits", "[particle]"
                   , char
                   , short, unsigned short
                   , int, unsigned int
                   , long, unsigned long
                   , long long, unsigned long long
                   ) {
  constexpr int Pos = sizeof(T) * 8 - 16;
  constexpr int Nbits = 16;
  TestType y = std::numeric_limits<TestType>::max();
  y &= ~( ~0uLL << Nbits );
  T s = 0;
  setbits<Pos,Nbits>( s, y );
  REQUIRE( s == ( static_cast<T>(y) << Pos ) );

  auto yy = getbits<Pos,Nbits,TestType>(s);
  REQUIRE( yy == y );
}

SCENARIO("setting species", "[particle]") {
  State s;
  REQUIRE( s.get<species>() == species::unknown );

  s.set(species::electron);
  REQUIRE( s.get<species>() == species::electron );
  s.set(species::positron);
  REQUIRE( s.get<species>() == species::positron );
  s.set(species::ion);
  REQUIRE( s.get<species>() == species::ion );
  s.set(species::photon);
  REQUIRE( s.get<species>() == species::photon );
}

SCENARIO("setting flags", "[particle]") {
  {
    State s;
    REQUIRE_FALSE( s.is(flag::exist));
  }

  auto f =
    []( flag X ) {
      State s;
      REQUIRE_FALSE(s.is(X));
      s.set(X);
      REQUIRE( s.is(X));
      s.reset(X);
      REQUIRE_FALSE( s.is(X));
    };

  for ( std::underlying_type_t<flag> i = 0; i < 8; ++i ) {
    f(static_cast<flag>(i));
  }
}

SCENARIO("setting birthplace", "[particle]") {
  State s;
  using bp_t = std::underlying_type_t<birthplace>;
  bp_t max = ( static_cast<bp_t>(1) << layout::size<birthplace>() );

  for ( bp_t i = 0; i < max; ++i ) {
    birthplace bp(i);
    s.set(bp);
    CAPTURE(i);
    CHECK(s.get<birthplace>() == i);
  }
}

// setting destination is in test_migration
