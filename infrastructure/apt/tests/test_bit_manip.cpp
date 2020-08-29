#include "apt/bit_manip.hpp"
#include "testfw/testfw.hpp"

using namespace apt;

TEMPLATE_TEST_CASE("Test setbits and getbits", "[apt]"
                   , char
                   , short, unsigned short
                   , int, unsigned int
                   , long, unsigned long
                   , long long, unsigned long long
                   ) {
  using T = long long;
  constexpr int Pos = sizeof(T) * 8 - 16;
  constexpr int Nbits = 16;
  TestType y = std::numeric_limits<TestType>::max();
  y &= ~( ~0uLL << Nbits );
  T s = 0;
  apt::setbits<Pos,Nbits>( s, y );
  REQUIRE( s == ( static_cast<T>(y) << Pos ) );

  auto yy = apt::getbits<Pos,Nbits>(s);
  REQUIRE( yy == y );
}
