#include "metric/log_spherical/geodesic.hpp"
#include "testfw/testfw.hpp"

using namespace apt;

#include "old_caldisp.h"

using std::cos;
using std::exp;
using std::log;
using std::sin;

auto PV = [](auto phi) {
  return phi - 2 * PI<double> * std::floor(phi / (2 * PI<double>));
};

SCENARIO("LogSpherical special cases", "[metric]") {
  aio::unif_real<double> unif;

  WHEN("velocity is radial") {
    int N = 100000;
    while (N--) {
      double dt = 0.001;
      Vec<double, 3> x_old(unif() * 3.0, unif() * (PI<double> - 0.2) + 0.1,
                           unif() * 2 * PI<double>);
      Vec<double, 3> p_old(unif() * 5.0, 0.0, 0.0);
      auto x = x_old;
      auto p = p_old;
      bool is_massive = unif() > 0.5;

      metric::LogSpherical::geodesic_move(x, p, dt, is_massive);

      CAPTURE(N, dt, x_old, p_old, x, p, is_massive);

      Vec<double, 3> v_old = p_old / std::sqrt(is_massive + apt::sqabs(p_old));
      REQUIRE(exp(x[0]) == Approx(v_old[0] * dt + exp(x_old[0])));
      REQUIRE(x[1] == Approx(x_old[1]));
      REQUIRE(x[2] == Approx(x_old[2]));

      for (int i = 0; i < 3; ++i)
        REQUIRE(p[i] == Approx(p_old[i]).margin(1e-12));
    }
  }

  WHEN("velocity is in meridional plane") {
    double dt = 0.001;
    int N = 100000;
    while (N--) {
      Vec<double, 3> x_old(unif() * 3.0, unif() * (PI<double> - 0.2) + 0.1,
                           unif() * 2 * PI<double>);
      Vec<double, 3> p_old(2.0 * (2 * unif() - 1), 2.0 * (2 * unif() - 1), 0.0);

      auto x = x_old;
      auto p = p_old;
      bool is_massive = unif() > 0.5;

      Coord::geodesic_move(x, p, dt, is_massive);

      CAPTURE(dt, x_old, p_old, x, p, is_massive);

      auto get_vzvx = [](double angle, double v_r, double v_t) {
        double v_z = cos(angle) * v_r - sin(angle) * v_t;
        double v_x = sin(angle) * v_r + cos(angle) * v_t;
        return std::array<double, 2>{v_z, v_x};
      };

      constexpr int R = 0;
      constexpr int THETA = 1;
      auto gamma = std::sqrt(is_massive + apt::sqabs(p_old));
      auto v_old = p_old / gamma;

      auto [v_z, v_x] = get_vzvx(x_old[THETA], v_old[R], v_old[THETA]);
      double z_new = exp(x_old[R]) * cos(x_old[THETA]) + v_z * dt;
      double x_new = exp(x_old[R]) * sin(x_old[THETA]) + v_x * dt;

      REQUIRE(exp(x[R]) * cos(x[THETA]) == Approx(z_new));
      REQUIRE(exp(x[R]) * sin(x[THETA]) == Approx(x_new));

      auto [p_z1, p_x1] = get_vzvx(x[THETA], p[R], p[THETA]);
      REQUIRE(p_z1 / gamma == Approx(v_z));
      REQUIRE(p_x1 / gamma == Approx(v_x));
    }
  }

  WHEN("location is at equator, velocity is in the plane") {
    double dt = 0.001;
    int N = 100000;
    while (N--) {
      Vec<double, 3> x_old(unif() * 3.0, PI<double> / 2.0,
                           unif() * 2 * PI<double>);
      Vec<double, 3> p_old(2.0 * (2 * unif() - 1), 0.0, 2.0 * (2 * unif() - 1));

      auto x = x_old;
      auto p = p_old;
      bool is_massive = unif() > 0.5;

      Coord::geodesic_move(x, p, dt, is_massive);

      CAPTURE(dt, x_old, p_old, x, p, is_massive);

      auto get_vxvy = [](double angle, double v_r, double v_p) {
        double v_x = cos(angle) * v_r - sin(angle) * v_p;
        double v_y = sin(angle) * v_r + cos(angle) * v_p;
        return std::array<double, 2>{v_x, v_y};
      };
      constexpr int R = 0;
      constexpr int PHI = 2;

      auto gamma = std::sqrt(is_massive + apt::sqabs(p_old));
      auto v_old = p_old / gamma;

      auto [v_x, v_y] = get_vxvy(x_old[PHI], v_old[R], v_old[PHI]);
      double x_new = exp(x_old[R]) * cos(x_old[PHI]) + v_x * dt;
      double y_new = exp(x_old[R]) * sin(x_old[PHI]) + v_y * dt;

      REQUIRE(exp(x[R]) * cos(x[PHI]) == Approx(x_new));
      REQUIRE(exp(x[R]) * sin(x[PHI]) == Approx(y_new));

      auto [p_x1, p_y1] = get_vxvy(x[PHI], p[R], p[PHI]);
      REQUIRE(p_x1 / gamma == Approx(v_x));
      REQUIRE(p_y1 / gamma == Approx(v_y));
    }
  }
}

SCENARIO(
    "LogSpherical, test against with old geodesic mover, which presumably is "
    "correct",
    "[metric]") {
  aio::unif_real<double> unif;

  int N = 1000000;
  const double dt = 0.001;
  while (N--) {
    // --- set up ----
    Vec<double, 3> x1(unif() * 3.0, unif() * PI<double>,
                      unif() * 2 * PI<double>);
    Vec<double, 3> p1(10.0 * (2 * unif() - 1), 10.0 * (2 * unif() - 1),
                      10.0 * (2 * unif() - 1));
    bool is_massive = unif() > 0.5;

    auto x2 = x1;
    Vec<double, 3> v2 = p1 / std::sqrt(is_massive + apt::sqabs(p1));

    const auto x_old = x1;
    // det is a determinant that plays a role in getting correct dx when axis
    // crossing happens
    auto det = (std::exp(x2[0]) + v2[0] * dt) * std::sin(x2[1]) +
               v2[1] * dt * std::cos(x2[1]);

    // --- run ----
    auto dx = metric::LogSpherical::geodesic_move(x1, p1, dt, is_massive);
    OLD_geodesic_move(x2, v2, dt);

    // --- test ---
    CAPTURE(dt, x1, p1, is_massive);

    REQUIRE(x1[0] == Approx(x2[0]));
    REQUIRE(x1[1] == Approx(x2[1]));

    REQUIRE(PV(x1[2]) == Approx(PV(x2[2])));

    Vec<double, 3> v1 = p1 / std::sqrt(is_massive + apt::sqabs(p1));
    REQUIRE(v1[0] == Approx(v2[0]));
    REQUIRE(v1[1] == Approx(v2[1]));
    REQUIRE(v1[2] == Approx(v2[2]));

    // test dx
    REQUIRE(dx[0] == Approx(x1[0] - x_old[0]));
    if (det > 0) {
      REQUIRE(dx[1] == Approx(x1[1] - x_old[1]));
      REQUIRE(dx[2] == Approx(x1[2] - x_old[2]));
    }
  }
}

SCENARIO("LogSpherical, test dx under axis crossing", "[metric]") {
  using Coord = mani::LogSphericalCoordSys;
  aio::unif_real<double> unif;

  int N = 100000;
  const double dt = 0.001;
  while (N) {
    // --- set up ----
    auto theta = unif() * dt * 0.1;
    if (unif() > 0.5) theta = PI<double> - theta;

    // NOTE r and pr are chosen to increase axis crossing
    Vec<double, 3> x1(unif() * 0.05, theta, unif() * 2 * PI<double>);
    Vec<double, 3> p1(1.0 * (2 * unif() - 1), 10.0 * (2 * unif() - 1),
                      10.0 * (2 * unif() - 1));
    bool is_massive = unif() > 0.5;

    const auto x_old = x1;
    Vec<double, 3> v_old = p1 / std::sqrt(is_massive + apt::sqabs(p1));
    // det is a determinant that plays a role in getting correct dx when axis
    // crossing happens
    auto det = (std::exp(x_old[0]) + v_old[0] * dt) * std::sin(x_old[1]) +
               v_old[1] * dt * std::cos(x_old[1]);
    if (det > 0)
      continue;
    else
      --N;

    // --- run ----
    auto dx = Coord::geodesic_move(x1, p1, dt, is_massive);

    // --- test ---
    CAPTURE(dt, x1, p1, is_massive);

    // test dx
    REQUIRE(dx[0] == Approx(x1[0] - x_old[0]));
    if (det > 0)
      REQUIRE(dx[1] == Approx(x1[1] - x_old[1]));
    else if (theta < PI<double> / 2.0)
      REQUIRE(dx[1] == Approx(-x1[1] - x_old[1]));
    else
      REQUIRE(dx[1] == Approx(2 * PI<double> - x1[1] - x_old[1]));

    REQUIRE(dx[2] == Approx(std::atan(v_old[2] * dt / det)));
  }
}
