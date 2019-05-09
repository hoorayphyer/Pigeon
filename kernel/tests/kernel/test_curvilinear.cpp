#include "all_in_one.hpp"
#include "kernel/curvilinear.hpp"

using namespace apt;

// SCENARIO("Cartesian", "[metric]") {
//   using Coord = knl::coord<knl::coordsys::Cartesian>;
//   util::Rng<double> rng;
//   rng.set_seed(std::time(0));

//   Vec<double,3> x_old( rng.uniform(-10,10), rng.uniform(-10,10), rng.uniform(-10,10) );
//   auto x = x_old;
//   Vec<double,3> v( rng.uniform(-10,10), rng.uniform(-10,10), rng.uniform(-10,10) );
//   double dt = rng.uniform();

//   Coord::geodesic_move( x, v, dt );
//   for ( int i = 0; i < 3; ++i )
//     REQUIRE( x[i] == Approx( v[i] * dt + x_old[i] ));
// }

#include "old_caldisp.h"

using std::sin;
using std::cos;
using std::exp;
using std::log;

auto PV = []( auto phi ) {
            return phi - 2 * PI<double> * std::floor( phi / ( 2 * PI<double>) );
          };

SCENARIO("LogSpherical special cases", "[metric]") {
  using Coord = knl::coord<knl::coordsys::LogSpherical>;
  aio::unif_real<double> unif;

  WHEN("velocity is radial") {
    int N = 10000;
    while(N--) {
      double dt = 0.001;
      Vec<double,3> x_old( unif() * 3.0, unif() * (PI<double>-0.2) + 0.1, unif() * 2 * PI<double> );
      Vec<double,3> p_old( unif() * 5.0, 0.0, 0.0 );
      auto x = x_old;
      auto p = p_old;
      bool is_massive = unif() > 0.5;

      Coord::geodesic_move( x, p, dt, is_massive );

      CAPTURE( dt, x_old, p_old, x, p, is_massive );

      auto v_old = p_old / std::sqrt( is_massive + apt::sqabs(p_old) );
      REQUIRE( exp(x[0]) == Approx( v_old[0] * dt + exp(x_old[0]) ));
      REQUIRE( x[1] == Approx( x_old[1] ) );
      REQUIRE( x[2] == Approx( x_old[2] ) );

      for ( int i = 0; i < 3; ++i )
        REQUIRE( p[i] == Approx(p_old[i]).margin(1e-12) );
    }
  }

  WHEN("velocity is in meridional plane") {
    double dt = 0.001;
    int N = 10000;
    while(N--) {
      Vec<double,3> x_old( unif() * 3.0, unif() * (PI<double>-0.2) + 0.1, unif() * 2 * PI<double> );
      Vec<double,3> p_old( 2.0 * ( 2 * unif() - 1 ), 2.0 * ( 2 * unif() - 1 ), 0.0 );

      auto x = x_old;
      auto p = p_old;
      bool is_massive = unif() > 0.5;

      Coord::geodesic_move( x, p, dt, is_massive );

      CAPTURE( dt, x_old, p_old, x, p, is_massive );

      auto get_vzvx =
        []( double angle, double v_r, double v_t ) {
          double v_z = cos(angle) * v_r - sin(angle) * v_t;
          double v_x = sin(angle) * v_r + cos(angle) * v_t;
          return std::array<double,2>{v_z, v_x};
        };

      constexpr int R = 0;
      constexpr int THETA = 1;
      auto gamma = std::sqrt( is_massive + apt::sqabs(p_old) );
      auto v_old = p_old / gamma;

      auto[v_z, v_x] = get_vzvx( x_old[THETA], v_old[R], v_old[THETA] );
      double z_new = exp(x_old[R]) * cos(x_old[THETA]) + v_z * dt;
      double x_new = exp(x_old[R]) * sin(x_old[THETA]) + v_x * dt;

      REQUIRE( exp(x[R]) * cos(x[THETA]) == Approx(z_new) );
      REQUIRE( exp(x[R]) * sin(x[THETA]) == Approx(x_new) );

      auto[p_z1, p_x1] = get_vzvx( x[THETA], p[R], p[THETA] );
      REQUIRE( p_z1 / gamma == Approx(v_z) );
      REQUIRE( p_x1 / gamma == Approx(v_x) );
    }

  }

  WHEN("location is at equator, velocity is in the plane") {
    double dt = 0.001;
    int N = 10000;
    while(N--) {
      Vec<double,3> x_old( unif() * 3.0, PI<double> / 2.0, unif() * 2 * PI<double> );
      Vec<double,3> p_old( 2.0 * ( 2 * unif() - 1 ), 0.0, 2.0 * ( 2 * unif() - 1 ) );

      auto x = x_old;
      auto p = p_old;
      bool is_massive = unif() > 0.5;

      Coord::geodesic_move( x, p, dt, is_massive );

      CAPTURE( dt, x_old, p_old, x, p, is_massive );


      auto get_vxvy =
        []( double angle, double v_r, double v_p ) {
          double v_x = cos(angle) * v_r - sin(angle) * v_p;
          double v_y = sin(angle) * v_r + cos(angle) * v_p;
          return std::array<double,2>{v_x, v_y};
        };
      constexpr int R = 0;
      constexpr int PHI = 2;

      auto gamma = std::sqrt( is_massive + apt::sqabs(p_old) );
      auto v_old = p_old / gamma;

      auto[v_x, v_y] = get_vxvy( x_old[PHI], v_old[R], v_old[PHI] );
      double x_new = exp(x_old[R]) * cos(x_old[PHI]) + v_x * dt;
      double y_new = exp(x_old[R]) * sin(x_old[PHI]) + v_y * dt;

      REQUIRE( exp(x[R]) * cos(x[PHI]) == Approx(x_new) );
      REQUIRE( exp(x[R]) * sin(x[PHI]) == Approx(y_new) );

      auto[p_x1, p_y1] = get_vxvy( x[PHI], p[R], p[PHI] );
      REQUIRE( p_x1 / gamma == Approx(v_x) );
      REQUIRE( p_y1 / gamma == Approx(v_y) );
    }
  }
}

SCENARIO("LogSpherical, test against with old geodesic mover, which presumably is correct", "[metric]") {
  using Coord = knl::coord<knl::coordsys::LogSpherical>;
  aio::unif_real<double> unif;

  int N = 100000;
  const double dt = 0.001;
  while (N--) {
    Vec<double,3> x1( unif() * 3.0, unif() * PI<double>, unif() * 2 * PI<double> );
    Vec<double,3> p1( 10.0 * ( 2 * unif() - 1 ), 10.0 * ( 2 * unif() - 1 ), 10.0 * ( 2 * unif() - 1 ) );
    bool is_massive = unif() > 0.5;

    auto x2 = x1;
    Vec<double,3> v2 = p1 / std::sqrt( is_massive + apt::sqabs(p1) );

    Coord::geodesic_move( x1, p1, dt, is_massive );
    OLD_geodesic_move( x2, v2, dt );

    CAPTURE( dt, x1, p1, is_massive );


    REQUIRE( x1[0] == Approx(x2[0]) );
    REQUIRE( x1[1] == Approx(x2[1]) );

    REQUIRE( PV(x1[2]) == Approx(PV(x2[2])));


    Vec<double,3> v1 = p1 / std::sqrt( is_massive + apt::sqabs(p1) );
    REQUIRE( v1[0] == Approx(v2[0]));
    REQUIRE( v1[1] == Approx(v2[1]));
    REQUIRE( v1[2] == Approx(v2[2]));

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
//     rdn[6*i + 1] = rng.uniform(0,PI<double>);
//     rdn[6*i + 2] = rng.uniform(0,2*PI<double>);

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

//       Coord::geodesic_move( x, v, dt );
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

//       logsph_geodesic_move_old( x, v, dt );
//     }
//     auto tf = high_resolution_clock::now();
//     auto dur = duration_cast<milliseconds>(tf - ti);
//     std::cout << "T2 = " << dur.count() << " ms" << std::endl;
//   }
// }
