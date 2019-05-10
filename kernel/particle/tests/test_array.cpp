#include "testfw/testfw.hpp"
#include "particle/array_impl.hpp"
#include "particle/particle.hpp"
#include "particle/virtual_particle.hpp"
#include "particle/c_particle.hpp"

using namespace particle;
using ptc_array = array<double,aio::Specs>;
using Ptc = Particle<double,aio::Specs>;
using vPtc = vParticle<double,aio::Specs>;
using cPtc = cParticle<double,aio::Specs>;
using Vec = apt::Vec<double,aio::Specs<double>::Dim>;

SCENARIO("array push_back and iterator", "[particle]") {
  ptc_array arr;

  WHEN("push back real particle") {
    Ptc ptc0 ( Vec(15,6,73), Vec(20,-3,-5), species::ion, flag::secondary );
    arr.push_back(ptc0);
    REQUIRE( arr.size() == 1 );
    AND_WHEN("by copy") {
      auto ptc = arr[0]; // NOTE ptc is a proxy hence no ref
      REQUIRE( ptc.q()[0] == 15 );
      REQUIRE( ptc.q()[1] == 6 );
      REQUIRE( ptc.q()[2] == 73 );

      REQUIRE( ptc.p()[0] == 20 );
      REQUIRE( ptc.p()[1] == -3 );
      REQUIRE( ptc.p()[2] == -5 );

      REQUIRE( ptc.is(species::ion) );
      REQUIRE( ptc.is(flag::secondary) );
    }

    Ptc ptc1 ( Vec(-45,10,-76), Vec(-13,-58,97), species::photon, flag::traced );
    arr.push_back(std::move(ptc1));
    REQUIRE( arr.size() == 2 );
    AND_WHEN("by move") {
      auto ptc = arr[1]; // NOTE ptc is a proxy hence no ref
      REQUIRE( ptc.q()[0] == -45 );
      REQUIRE( ptc.q()[1] == 10 );
      REQUIRE( ptc.q()[2] == -76 );

      REQUIRE( ptc.p()[0] == -13 );
      REQUIRE( ptc.p()[1] == -58 );
      REQUIRE( ptc.p()[2] == 97 );

      REQUIRE( ptc.is(species::photon) );
      REQUIRE( ptc.is(flag::traced) );
    }
  }

  WHEN("push back virtual particle") {
    apt::array<double,3> q { 15,6,73 };
    apt::array<double,3> p { 20,-3,-5 };
    unsigned long long state = 147;
    vPtc ptc0 ( q, p, state );
    arr.push_back(std::move(ptc0));
    REQUIRE( arr.size() == 1 );
    // NOTE no push_back by copy for virtual particle in order to preserve semantics
    AND_WHEN("by move") {
      auto ptc = arr[0]; // NOTE ptc is a proxy hence no ref
      REQUIRE( ptc.q()[0] == 15 );
      REQUIRE( ptc.q()[1] == 6 );
      REQUIRE( ptc.q()[2] == 73 );

      REQUIRE( ptc.p()[0] == 20 );
      REQUIRE( ptc.p()[1] == -3 );
      REQUIRE( ptc.p()[2] == -5 );

      REQUIRE( ptc.state() == 147 );
    }
  }

  WHEN("push back c particle") {
    cPtc ptc0 ( Ptc(Vec(15,6,73), Vec(20,-3,-5), species::ion, flag::secondary) );
    arr.push_back(ptc0);
    REQUIRE( arr.size() == 1 );
    AND_WHEN("by copy") {
      auto ptc = arr[0]; // NOTE ptc is a proxy hence no ref
      REQUIRE( ptc.q()[0] == 15 );
      REQUIRE( ptc.q()[1] == 6 );
      REQUIRE( ptc.q()[2] == 73 );

      REQUIRE( ptc.p()[0] == 20 );
      REQUIRE( ptc.p()[1] == -3 );
      REQUIRE( ptc.p()[2] == -5 );

      REQUIRE( ptc.is(species::ion) );
      REQUIRE( ptc.is(flag::secondary) );
    }

    cPtc ptc1 ( Ptc(Vec(-45,10,-76), Vec(-13,-58,97), species::photon, flag::traced) );
    arr.push_back(std::move(ptc1));
    REQUIRE( arr.size() == 2 );
    AND_WHEN("by move") {
      auto ptc = arr[1]; // NOTE ptc is a proxy hence no ref
      REQUIRE( ptc.q()[0] == -45 );
      REQUIRE( ptc.q()[1] == 10 );
      REQUIRE( ptc.q()[2] == -76 );

      REQUIRE( ptc.p()[0] == -13 );
      REQUIRE( ptc.p()[1] == -58 );
      REQUIRE( ptc.p()[2] == 97 );

      REQUIRE( ptc.is(species::photon) );
      REQUIRE( ptc.is(flag::traced) );
    }
  }
}

SCENARIO("array back inserter", "[particle]") {
  ptc_array arr;
  auto inserter = std::back_inserter(arr);
  (*inserter++) = Ptc( Vec(15,6,73), Vec(20,-3,-5), species::ion, flag::secondary );
  REQUIRE( arr.size() == 1 );

  auto ptc = arr[0]; // NOTE ptc is a proxy hence no ref
  REQUIRE( ptc.q()[0] == 15 );
  REQUIRE( ptc.q()[1] == 6 );
  REQUIRE( ptc.q()[2] == 73 );

  REQUIRE( ptc.p()[0] == 20 );
  REQUIRE( ptc.p()[1] == -3 );
  REQUIRE( ptc.p()[2] == -5 );

  REQUIRE( ptc.is(species::ion) );
  REQUIRE( ptc.is(flag::secondary) );
}

SCENARIO("array erase", "[particle]") {
  ptc_array arr;
  Ptc ptc0 ( Vec(15,6,73), Vec(20,-3,-5), species::ion, flag::secondary );
  for ( int i = 0; i < 50; ++i )
    arr.push_back(ptc0);
  THEN("erase will mark particles as empty but array size is unchanged yet") {
    WHEN("erase range is normal") {
      for ( int i = 10; i < 20; ++i )
        REQUIRE_FALSE( arr[i].is(flag::empty) );

      arr.erase( 10, 20 );

      for ( int i = 10; i < 20; ++i )
        REQUIRE( arr[i].is(flag::empty) );
      REQUIRE( arr.size() == 50 );
    }

    WHEN("erase range is inverted") {
      for ( int i = 40; i > 30; --i )
        REQUIRE_FALSE( arr[i].is(flag::empty) );

      arr.erase( 40, 30 );

      for ( int i = 40; i > 30; --i )
        REQUIRE( arr[i].is(flag::empty) );
      REQUIRE( arr.size() == 50 );
    }

    AND_THEN("automatic off-bound protection") {
      arr.erase( 40, 6000 );
      REQUIRE( arr.size() == 50 );
    }
  }
}

SCENARIO("array resize", "[particle]") {
  ptc_array arr;
  Ptc ptc0 ( Vec(15,6,73), Vec(20,-3,-5), species::ion, flag::secondary );
  for ( int i = 0; i < 50; ++i )
    arr.push_back(ptc0);

  WHEN("resize to shrink") {
    arr.resize( 30 );
    REQUIRE(arr.size() == 30 );
  }

  WHEN("resize to expand") {
    arr.resize( 80 );
    REQUIRE(arr.size() == 80 );
    for ( int i = 0; i < 50; ++i )
      REQUIRE_FALSE( arr[i].is(flag::empty) );

    for ( int i = 50; i < 80; ++i )
      REQUIRE( arr[i].is(flag::empty) );
  }

}
