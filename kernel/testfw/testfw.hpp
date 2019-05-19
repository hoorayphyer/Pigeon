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

  template < typename T, template < typename > class Distribution >
  class Dist : public Distribution<T> {
  private:
    eng_t _eng {now()};

  public:
    using Distribution<T>::Distribution;

    inline void seed( std::size_t s ) { _eng.seed(s); }

    inline T operator() () {
      return static_cast<Distribution<T>&>(*this)(_eng);
    }
  };

  template < typename T = int >
  using unif_int = Dist<T,std::uniform_int_distribution>; // [lower,upper]

  template < typename T = double >
  using unif_real = Dist<T,std::uniform_real_distribution>;

  template < typename T = double >
  using gauss_real = Dist<T,std::normal_distribution>;
}

namespace mpi { struct CartComm; }

// define Specs so that particle struct can be complete
#include "apt/type_traits.hpp"
namespace aio {
  template < typename T >
  struct Specs {
    using value_type = T;
    static constexpr int Dim = 3;
    using state_type = apt::copy_cvref_t<T,unsigned long long>;

    static_assert( 8 * sizeof( state_type ) >= 64 );
  };
}

namespace aio {
  template < int... I>
  struct IndexType {
    constexpr static apt::array<int, sizeof...(I)> get() noexcept {return {I...};}
  };
}

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

  // NOTE Notation: XxYxZ is the cartesian partition. Each one of X,Y,Z can be positive ( meaning nonperiodic ) or negative ( meaning periodic ).
  template < int Size, typename CartComm = mpi::CartComm, typename WorldComm  >
  auto make_cart( apt::array<int,Size> topo, const WorldComm& world ) {
    std::vector<int> cart_dims;
    std::vector<bool> periodic;
    for ( auto i : topo ) {
      bool is_neg = (i < 0);
      cart_dims.push_back( is_neg ? -i : i );
      periodic.push_back( is_neg );
    }
    return make_cart( std::move(cart_dims), std::move(periodic), world );
  }
}

// #include "logger/logger.hpp"
// namespace aio {
  
// }

#endif
