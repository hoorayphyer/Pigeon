#include <iostream>
#include "kernel/coordinate.hpp"
#include "apt/vec.hpp"
#include "apt/print_vec.hpp"
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
  for ( int i = 0; i < 3; ++i )
    REQUIRE( dx[i] == Approx( v[i] * dt));

  for ( int i = 0; i < 3; ++i )
    REQUIRE( x[i] == Approx( v[i] * dt + x_old[i] ));
}

#include "old_caldisp.h"

using std::sin;
using std::cos;
using std::exp;
using std::log;

constexpr double PI = PI_CONST;

SCENARIO("LogSpherical", "[knl]") {
  using Coord = knl::coord<knl::coordsys::LogSpherical>;
  util::Rng<double> rng;
  rng.set_seed(std::time(0));

  GIVEN("NOT on axes") {
    WHEN("velocity is radial") {
      Vec<double,3> x_old( rng.uniform(0,3), rng.uniform(0.1,PI-0.1), rng.uniform(0,2*PI) );
      Vec<double,3> v_old( rng.uniform(0,5), 0.0, 0.0 );
      v_old /= std::sqrt( 1.0 + apt::sqabs(v_old) );
      auto x = x_old;
      auto v = v_old;
      double dt = rng.uniform( 0.0, 0.1 );

      Vec<double,3> dx = Coord::geodesic_move( x, v, dt );

      CAPTURE( dt, x_old, v_old, x, v, dx );

      REQUIRE( exp(x[0]) == Approx( v[0] * dt + exp(x_old[0]) ));
      REQUIRE( x[1] == Approx( x_old[1] ) );
      REQUIRE( x[2] == Approx( x_old[2] ) );

      for ( int i = 0; i < 3; ++i )
        REQUIRE( v[i] == Approx(v_old[i]).margin(1e-12) );

      REQUIRE( dx[0] == Approx( x[0] - x_old[0] ) );
      REQUIRE( dx[1] == Approx( 0.0 ).margin(1e-12) );
      REQUIRE( dx[2] == Approx( 0.0 ).margin(1e-12) );
    }

    WHEN("velocity is in meridional plane") {
      Vec<double,3> x_old( rng.uniform(0,3), rng.uniform(0.1,PI-0.1), 0.0 );
      Vec<double,3> v_old( rng.uniform(-2,2) * 0, rng.uniform(-2,2), 0.0 );
      v_old /= std::sqrt( 1.0 + apt::sqabs(v_old) );

      auto x = x_old;
      auto v = v_old;
      double dt = rng.uniform( 0.0, 0.1 );

      Vec<double,3> dx = Coord::geodesic_move( x, v, dt );

      CAPTURE( dt, x_old, v_old, x, v );

      auto get_vzvx =
        []( double angle, double v_r, double v_t ) {
          double v_z = cos(angle) * v_r - sin(angle) * v_t;
          double v_x = sin(angle) * v_r + cos(angle) * v_t;
          return std::array<double,2>{v_z, v_x};
        };

      constexpr int R = 0;
      constexpr int THETA = 1;
      auto[v_z, v_x] = get_vzvx( x_old[THETA], v_old[R], v_old[THETA] );
      double z_new = exp(x_old[R]) * cos(x_old[THETA]) + v_z * dt;
      double x_new = exp(x_old[R]) * sin(x_old[THETA]) + v_x * dt;

      REQUIRE( exp(x[R]) * cos(x[THETA]) == Approx(z_new) );
      REQUIRE( exp(x[R]) * sin(x[THETA]) == Approx(x_new) );

      auto[v_z1, v_x1] = get_vzvx( x[THETA], v[R], v[THETA] );
      REQUIRE( v_z1 == Approx(v_z) );
      REQUIRE( v_x1 == Approx(v_x) );

    }

    WHEN("location is at equator, velocity is in the plane") {
      Vec<double,3> x_old( rng.uniform(0,3), PI / 2.0, rng.uniform(0.1, 2*PI) );
      Vec<double,3> v_old( rng.uniform(-2,2), 0.0, rng.uniform(-2,2) );
      v_old /= std::sqrt( 1.0 + apt::sqabs(v_old) );

      auto x = x_old;
      auto v = v_old;
      double dt = rng.uniform( 0.0, 0.1 );

      Vec<double,3> dx = Coord::geodesic_move( x, v, dt );

      CAPTURE( dt, x_old, v_old, x, v );

      auto get_vxvy =
        []( double angle, double v_r, double v_p ) {
          double v_x = cos(angle) * v_r - sin(angle) * v_p;
          double v_y = sin(angle) * v_r + cos(angle) * v_p;
          return std::array<double,2>{v_x, v_y};
        };
      constexpr int R = 0;
      constexpr int PHI = 2;

      auto[v_x, v_y] = get_vxvy( x_old[PHI], v_old[R], v_old[PHI] );
      double x_new = exp(x_old[R]) * cos(x_old[PHI]) + v_x * dt;
      double y_new = exp(x_old[R]) * sin(x_old[PHI]) + v_y * dt;

      REQUIRE( exp(x[R]) * cos(x[PHI]) == Approx(x_new) );
      REQUIRE( exp(x[R]) * sin(x[PHI]) == Approx(y_new) );

      auto[v_x1, v_y1] = get_vxvy( x[PHI], v[R], v[PHI] );
      REQUIRE( v_x1 == Approx(v_x) );
      REQUIRE( v_y1 == Approx(v_y) );

    }
  }

  // TODO test on axes

  SECTION("compare with old geodesic mover, which presumably is correct") {
    util::Rng<double> rng;
    rng.set_seed(std::time(0));

    int N = 100;
    for ( int i = 0; i < N; ++i ) {
      const double dt = rng.uniform(0.0, 0.01);
      Vec<double,3> x1( rng.uniform(0,3), rng.uniform(0, PI), rng.uniform(0,2*PI) );
      Vec<double,3> v1( rng.uniform(-10,10), rng.uniform(-10,10), rng.uniform(-10,10) );
      v1 /= std::sqrt( 1.0 + apt::sqabs(v1) );

      auto x2 = x1;
      auto v2 = v1;

      CAPTURE(dt, x1, v1);

      Vec<double,3> dx1 = Coord::geodesic_move( x1, v1, dt );
      Vec<double,3> dx2 = OLD_geodesic_move( x2, v2, dt );

      CAPTURE(dx1, dx2);
      REQUIRE( x1[0] == Approx(x2[0]));
      REQUIRE( dx1[0] == Approx(dx2[0]));
      REQUIRE( x1[1] == Approx(x2[1]));
      REQUIRE( dx1[1] == Approx(dx2[1]));

      // REQUIRE( x1[2] == Approx(x2[2]));
      // REQUIRE( dx1[2] == Approx(dx2[2]));

      // REQUIRE( v1[0] == Approx(v2[0]));
      // REQUIRE( v1[1] == Approx(v2[1]));
      // REQUIRE( v1[2] == Approx(v2[2]));

    }

  }
}

#include <chrono>
using namespace std::chrono;

// SCENARIO("measure performance of LogSpherical", "[knl]") {
//   using Coord = knl::coord<knl::coordsys::LogSpherical>;
//   util::Rng<double> rng;
//   rng.set_seed(std::time(0));

//   const double dt = 0.001;
//   int MIL = 1000*1000;
//   std::vector<double> rdn ( MIL * 6 );
//   for ( int i = 0; i < MIL; ++i ) {
//     rdn[6*i + 0] = rng.uniform(0,3);
//     rdn[6*i + 1] = rng.uniform(0,PI);
//     rdn[6*i + 2] = rng.uniform(0,2*PI);

//     rdn[6*i + 3] = rng.uniform(-10,10);
//     rdn[6*i + 4] = rng.uniform(-10,10);
//     rdn[6*i + 5] = rng.uniform(-10,10);

//     double gamma
//       = std::sqrt( 1.0 + rdn[6*i + 3] * rdn[6*i + 3]
//                    + rdn[6*i + 4] * rdn[6*i + 4]
//                    + rdn[6*i + 5] * rdn[6*i + 5]  );

//     rdn[6*i + 3] /= gamma;
//     rdn[6*i + 4] /= gamma;
//     rdn[6*i + 5] /= gamma;
//   }

//   int X = 100;
//   { // new geodesic_move
//     auto ti = high_resolution_clock::now();
//     for ( int I = 0; I < X * MIL; ++I ) {
//       int i = I % MIL;
//       Vec<double,3> x( rdn[6*i + 0], rdn[6*i + 1], rdn[6*i + 2] );
//       Vec<double,3> v( rdn[6*i + 3], rdn[6*i + 4], rdn[6*i + 5] );

//       Vec<double,3> dx = Coord::geodesic_move( x, v, dt );
//     }
//     auto tf = high_resolution_clock::now();
//     auto dur = duration_cast<milliseconds>(tf - ti);
//     std::cout << "T1 = " << dur.count() << " ms" << std::endl;
//   }

//   { // the old caldisp
//     auto ti = high_resolution_clock::now();
//     for ( int I = 0; I < X * MIL; ++I ) {
//       int i = I % MIL;
//       Vec3<double> x( rdn[6*i + 0], rdn[6*i + 1], rdn[6*i + 2] );
//       Vec3<double> v( rdn[6*i + 3], rdn[6*i + 4], rdn[6*i + 5] );

//       Vec3<double> dx3 = logsph_geodesic_move_old( x, v, dt );
//     }
//     auto tf = high_resolution_clock::now();
//     auto dur = duration_cast<milliseconds>(tf - ti);
//     std::cout << "T2 = " << dur.count() << " ms" << std::endl;
//   }
// }
