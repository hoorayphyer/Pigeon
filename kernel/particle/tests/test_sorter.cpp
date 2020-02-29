#include "testfw/testfw.hpp"
#include "particle/sorter.hpp"
#include "particle/array_impl.hpp"

#include <random>
#include <algorithm>

using namespace particle;

template < typename T >
struct S {
  using value_type = T;
  static constexpr int Dim = 3;
  using state_type = apt::copy_cvref_t<T,unsigned long long>;
};

SCENARIO("Test erase nonexist", "[particle]") {
  using Real = double;
  array<Real,S> ptcs;
  WHEN("empty buffer") {
    erase_nonexist(ptcs);
    REQUIRE(ptcs.size() == 0);
  }
  WHEN("buffer has one ptc") {
    ptcs.emplace_back({},{},1.0,species::electron);
    REQUIRE(ptcs.size() == 1);
    REQUIRE(ptcs[0].is(flag::exist));
    erase_nonexist(ptcs);
    REQUIRE(ptcs.size() == 1);
    REQUIRE(ptcs[0].is(flag::exist));
  }
  WHEN("buffer has only empty particles") {
    for ( int i = 0; i < 100; ++i )
      ptcs.emplace_back();
    erase_nonexist(ptcs);
    REQUIRE(ptcs.size() == 0 );
  }
  WHEN("buffer has only nonempty particles") {
    for ( int i = 0; i < 100; ++i )
      ptcs.emplace_back({},{},1.0);
    erase_nonexist(ptcs);
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

      erase_nonexist(ptcs);
      REQUIRE(ptcs.size() == num_ptcs - num_empties);
      for ( const auto& x : ptcs )
        REQUIRE(x.is(flag::exist));
    }

  }
}

SCENARIO("Test permutation", "[particle]") {
  std::size_t N = 1000 * 1000;
  const auto data =
    [N]() {
      std::vector<std::size_t> data(N);
      for ( std::size_t i = 0; i < data.size(); ++i ) data[i] = i;
      return data;
    }();

  std::random_device rd;
  std::mt19937 g(rd());

  int M = 100;
  while (M--) {
    auto d = data;
    auto f = [&d]( std::size_t i, std::size_t j ) mutable {
               std::swap(d[i],d[j]);
             };

    auto perm = data;
    std::shuffle(perm.begin(), perm.end(), g);
    const auto perm_copy = perm;

    permute(perm,f);
    for ( std::size_t i = 0; i < d.size(); ++i ) {
      CHECK( d[i] == perm_copy[i] );
      CHECK( perm[i] == i );
    }

  };
}
