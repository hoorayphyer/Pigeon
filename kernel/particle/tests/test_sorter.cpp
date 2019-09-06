#include "testfw/testfw.hpp"
#include "particle/sorter.hpp"
#include "particle/array_impl.hpp"

using namespace particle;

template < typename T >
struct S {
  using value_type = T;
  static constexpr int Dim = 3;
  using state_type = apt::copy_cvref_t<T,unsigned long long>;
};

SCENARIO("Test particle sorter", "[particle]") {
  using Real = double;
  array<Real,S> ptcs;
  WHEN("empty buffer") {
    sort(ptcs);
    REQUIRE(ptcs.size() == 0);
  }
  WHEN("buffer has one ptc") {
    ptcs.emplace_back({},{},1.0,species::electron);
    REQUIRE(ptcs.size() == 1);
    REQUIRE(ptcs[0].is(flag::exist));
    sort(ptcs);
    REQUIRE(ptcs.size() == 1);
    REQUIRE(ptcs[0].is(flag::exist));
  }
  WHEN("buffer has only empty particles") {
    for ( int i = 0; i < 100; ++i )
      ptcs.emplace_back();
    sort(ptcs);
    REQUIRE(ptcs.size() == 0 );
  }
  WHEN("buffer has only nonempty particles") {
    for ( int i = 0; i < 100; ++i )
      ptcs.emplace_back({},{},1.0);
    sort(ptcs);
    REQUIRE(ptcs.size() == 100 );
  }
  WHEN("random test") {
    aio::unif_int<int> unif(1, 1000);
    int N = 1000;
    while ( N-- ) {
      ptcs.resize(0);
      const int num_ptcs = unif();
      int num_empties = 0;
      for ( int i = 0; i < num_ptcs; ++i ) {
        bool is_empty = (unif() % 2 == 0);
        if ( is_empty ) {
          ++num_empties;
          ptcs.emplace_back();
        } else {
          ptcs.emplace_back({},{},1.0,species::electron);
        }
      }

      sort(ptcs);
      REQUIRE(ptcs.size() == num_ptcs - num_empties);
      for ( const auto& x : ptcs )
        REQUIRE(x.is(flag::exist));
    }

  }
}
