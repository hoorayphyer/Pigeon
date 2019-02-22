#include <iostream>
#include "kernel/coordinate.hpp"
#include "apt/vec.hpp"
#include "catch2/catch.hpp"
#include "utility/rng.hpp"
#include <ctime>


using namespace apt;

SCENARIO("Cartesian", "[knl]") {
  using Coord = knl::coord<knl::coordsys::Cartesian>;
  util::Rng<double> rng;
  rng.set_seed(std::time(0));

  Vec<double,3> x_old( rng.uniform(-10,10), rng.uniform(-10,10), rng.uniform(-10,10) );
  auto x = x_old;
  Vec<double,3> v( rng.uniform(-10,10), rng.uniform(-10,10), rng.uniform(-10,10) );
  double dt = rng.uniform();

  Vec<double,3> dx = Coord::geodesic_move( x, v, dt );
  REQUIRE( std::get<0>(dx) == Approx(std::get<0>(v) * dt));
  REQUIRE( std::get<1>(dx) == Approx(std::get<1>(v) * dt));
  REQUIRE( std::get<2>(dx) == Approx(std::get<2>(v) * dt));

  REQUIRE( std::get<0>(x) == Approx(std::get<0>(v) * dt + std::get<0>(x_old)));
  REQUIRE( std::get<1>(x) == Approx(std::get<1>(v) * dt + std::get<1>(x_old)));
  REQUIRE( std::get<2>(x) == Approx(std::get<2>(v) * dt + std::get<2>(x_old)));
}

// TODO do more tests
SCENARIO("LogSpherical", "[knl]") {
  using Coord = knl::coord<knl::coordsys::LogSpherical>;
  util::Rng<double> rng;
  rng.set_seed(std::time(0));

  GIVEN("velocity is radial") {
    WHEN("NOT on axes") {
      Vec<double,3> x_old( rng.uniform(0,3), rng.uniform(0.1,PI<double>-0.1), rng.uniform(0,2*PI<double>) );
      Vec<double,3> v_old( rng.uniform(0,5), 0.0, 0.0 );
      v_old /= std::sqrt( 1.0 + apt::sqabs(v_old) );
      auto x = x_old;
      auto v = v_old;
      double dt = rng.uniform( 0.0, 0.1 );

      Vec<double,3> dx = Coord::geodesic_move( x, v, dt );
      REQUIRE( std::exp(std::get<0>(x)) == Approx(std::get<0>(v) * dt + std::exp(std::get<0>(x_old)) ));
      REQUIRE( std::get<1>(x) == Approx(std::get<1>(x_old)));
      REQUIRE( std::get<2>(x) == Approx(std::get<2>(x_old)));

      REQUIRE( std::get<0>(v) == Approx(std::get<0>(v_old)));
      // TODO for v_theta, comparing with 0.0 gives artificial fail
      // std::cout << " v = " << std::get<1>(v) << ", v_old = " << std::get<1>(v_old) << std::endl;
      // REQUIRE( std::get<1>(v) == Approx(std::get<1>(v_old)).margin(1e-12));
      REQUIRE( std::get<2>(v) == Approx(std::get<2>(v_old)));

      REQUIRE( std::get<0>(dx) == Approx(std::get<0>(x) - std::get<0>(x_old)) );
      // REQUIRE( std::get<1>(dx) == Approx(0.0) );
      // REQUIRE( std::get<2>(dx) == Approx(0.0) );
    }

  }

  // SECTION("compare with old geodesic mover, which presumably is correct") {
  //   util::Rng<double> rng;
  //   rng.set_seed(std::time(0));

  //   int N = 10;
  //   for ( int i = 0; i < N; ++i ) {
  //     const double dt = rng.uniform(0.0, 0.01);
  //     Vec<double,3> x1( rng.uniform(0,3), rng.uniform(0, PI<double>), rng.uniform(0,2*PI<double>) );
  //     Vec<double,3> v1( rng.uniform(-10,10), rng.uniform(-10,10), rng.uniform(-10,10) );
  //     v1 /= std::sqrt( 1.0 + apt::sqabs(v1) );

  //     auto x2 = x1;
  //     auto v2 = v1;

  //     Vec<double,3> dx1 = Coord::geodesic_move( x1, v1, dt );
  //     Vec<double,3> dx2 = logsph_geodesic_move_old( x2, v2, dt );

  //     REQUIRE( std::get<0>(x1) == Approx(std::get<0>(x2)));
  //     REQUIRE( std::get<1>(x1) == Approx(std::get<1>(x2)));
  //     REQUIRE( std::get<2>(x1) == Approx(std::get<2>(x2)));

  //     REQUIRE( std::get<0>(v1) == Approx(std::get<0>(v2)));
  //     REQUIRE( std::get<1>(v1) == Approx(std::get<1>(v2)));
  //     REQUIRE( std::get<2>(v1) == Approx(std::get<2>(v2)));

  //     REQUIRE( std::get<0>(dx1) == Approx(std::get<0>(dx2)));
  //     REQUIRE( std::get<1>(dx1) == Approx(std::get<1>(dx2)));
  //     REQUIRE( std::get<2>(dx1) == Approx(std::get<2>(dx2)));
  //   }

  // }
}

#include "old_caldisp.h"
SCENARIO("measure performance of LogSpherical", "[knl]") {
  using Coord = knl::coord<knl::coordsys::LogSpherical>;
  util::Rng<double> rng;
  rng.set_seed(std::time(0));

  const double dt = 0.001;
  int N = 1000*1000;
  std::vector<double> rdn ( N * 6 );
  for ( int i = 0; i < N; ++i ) {
    rdn[6*i + 0] = rng.uniform(0,3);
    rdn[6*i + 1] = rng.uniform(0,PI<double>);
    rdn[6*i + 2] = rng.uniform(0,2*PI<double>);

    rdn[6*i + 3] = rng.uniform(-10,10);
    rdn[6*i + 4] = rng.uniform(-10,10);
    rdn[6*i + 5] = rng.uniform(-10,10);
  }

  int M = 1000;

  { // new geodesic_move
    std::time_t ti = std::time(nullptr);
    for ( int I = 0; I < M * N; ++I ) {
      int i = I % N;
      Vec<double,3> x1( rdn[6*i + 0], rdn[6*i + 1], rdn[6*i + 2] );
      Vec<double,3> v1( rdn[6*i + 3], rdn[6*i + 4], rdn[6*i + 5] );
      v1 /= std::sqrt( 1.0 + apt::sqabs(v1) );

      Vec<double,3> dx1 = Coord::geodesic_move( x1, v1, dt );
    }
    std::time_t tf = std::time(nullptr);
    std::cout << "T1 = " << tf - ti << " s" << std::endl;
  }

  { //old caldisp optimized
    std::time_t ti = std::time(nullptr);
    for ( int I = 0; I < 1000 * N; ++I ) {
      int i = I % N;
      Vec<double,3> x2( rdn[6*i + 0], rdn[6*i + 1], rdn[6*i + 2] );
      Vec<double,3> v2( rdn[6*i + 3], rdn[6*i + 4], rdn[6*i + 5] );
      v2 /= std::sqrt( 1.0 + apt::sqabs(v2) );

      Vec<double,3> dx2 = logsph_geodesic_move_old_modified( x2, v2, dt );
    }
    std::time_t tf = std::time(nullptr);
    std::cout << "T2 = " << tf - ti << " s" << std::endl;
  }

  { // the old caldisp
    std::time_t ti = std::time(nullptr);
    for ( int I = 0; I < 1000 * N; ++I ) {
      int i = I % N;
      Vec3<double> x3( rdn[6*i + 0], rdn[6*i + 1], rdn[6*i + 2] );
      Vec3<double> v3( rdn[6*i + 3], rdn[6*i + 4], rdn[6*i + 5] );
      v3 /= std::sqrt( 1.0 + v3.dot(v3) );

      Vec3<double> dx3 = logsph_geodesic_move_old( x3, v3, dt );
    }
    std::time_t tf = std::time(nullptr);
    std::cout << "T3 = " << tf - ti << " s" << std::endl;
  }
}
