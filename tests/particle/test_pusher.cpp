#include "particle/pusher.cpp"
#include "old_vay_pusher.h"
#include "utility/rng.hpp"
#include "catch2/catch.hpp"
#include <ctime>

using namespace particle;

SCENARIO("Test against old Vay Pusher") {
  util::Rng<double> rng;
  rng.set_seed(std::time(0));
  const int N = 10000000;
  for ( int i = 0; i < N; ++i ) {
    apt::Vec<double, 3> E( rng.uniform(-1000, 1000), rng.uniform(-1000, 1000), rng.uniform(-1000, 1000) );
    apt::Vec<double, 3> B( rng.uniform(-1000, 1000), rng.uniform(-1000, 1000), rng.uniform(-1000, 1000) );
    apt::Vec<double, 3> p( rng.uniform(-10, 10), rng.uniform(-10, 10), rng.uniform(-10, 10) );

    double lambda = 0.01;// lambda = charge_to_mass * dt;

    auto dp_old = lorentz_force_old( lambda, p, E, B );
    auto dp_new = force::lorentz( lambda, p, E, B );
    for ( int i = 0; i < 3; ++i ) {
      REQUIRE( dp_old[i] == Approx(dp_new[i]) );
    }
  }

}
