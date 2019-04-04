#ifndef _TEST_ALL_IN_ONE_HPP_
#define _TEST_ALL_IN_ONE_HPP_

#include "apt/print.hpp"
#include "catch2/catch.hpp"
#include <random>
#include <ctime>
#include <cmath>
#include <iostream>
namespace aio {
  using eng_t = std::default_random_engine;
  auto now() { return std::time(0); }

  template < typename T = int >
  using unif_int = std::uniform_int_distribution<T>; // [lower,upper]

  template < typename T = double >
  using unif_real = std::uniform_real_distribution<T>;
}

namespace mpi { struct CartComm; }

namespace aio {

  template < typename CartComm = mpi::CartComm, typename WorldComm  >
  auto make_cart( std::vector<int> dims, std::vector<bool> periodic, const WorldComm& world ) {
    std::optional<CartComm> cart;
    int size = 1;
    for ( auto x : dims ) size *= x;
    if ( world.size() >= size ) {
      auto comm = world.split( (world.rank() < size) );
      if ( world.rank() < size )
        cart.emplace( *comm, dims, periodic );
    }

    return cart;
  }
}

#include <array>
namespace aio {
  template < int... I>
  struct IndexType {
    constexpr static std::array<int, sizeof...(I)> get() noexcept {return {I...};}
  };
}

#endif
