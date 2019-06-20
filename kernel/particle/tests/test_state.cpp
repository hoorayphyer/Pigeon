#include "testfw/testfw.hpp"
#include "particle/state.hpp"

using namespace particle;

using T = unsigned long long;

struct State : public StateExpression<State,T> {
private:
  T _state = 0;
public:
  constexpr T& state() noexcept { return _state; }
  constexpr T state() const noexcept { return _state; }
};

SCENARIO("setting species", "[particle]") {
  State s;

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
  State s;
  s.set(flag::empty);
  REQUIRE( s.is(flag::empty));
  s.set(flag::traced);
  REQUIRE( s.is(flag::traced));
}

// TODO fix tracing, interface too bad
SCENARIO("setting tracing", "[particle]") {
  State s;
  birthplace bp(265);
  // serial_number sn(5654654);
  s.set(bp);
  // s.set(sn);
  REQUIRE(s.get<birthplace>() == 265);
  // REQUIRE(s.get<serial_number>() == 5654654);
}

