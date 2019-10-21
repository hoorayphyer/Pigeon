#include "testfw/testfw.hpp"
#include "manifold/curvilinear.hpp"

// SCENARIO("measure performance of LogSpherical", "[mani]") {
//   using Coord = mani::coord<mani::coordsys::LogSpherical>;
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
