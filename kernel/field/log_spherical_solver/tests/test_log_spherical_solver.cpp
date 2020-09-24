#include "testfw/testfw.hpp"
#include "field/log_spherical_solver/updater_impl.hpp"

using namespace field;
using R = double; // NOTE the following error bound doesn't work for float
constexpr int D = 2;

TEST_CASE("test curlE", "[field]") {
  const int N = 128;
  const apt::Grid<R,D> grid = {{{0.0, std::log(30.0), N}, {0.0, std::acos(-1.0), N}}};
  constexpr int g = 1;

  auto range = apt::make_range<D>({}, {grid[0].dim(), grid[1].dim()}, g);

  std::optional<int> surf(0);
  std::optional<int> outer(N);
  std::optional<int> lo_axis(0);
  std::optional<int> hi_axis(N+1);

  const int b[D] = {0, 0};
  const int e[D] = {N, N+1};

  Diff_in_CurlE<R, D> diff_r{0, surf, outer, grid[r_].delta()};
  Diff_in_CurlE<R, D> diff_th{1, surf, outer, grid[th_].delta()};

  Field<R, 3, D> E = {range};
  Field<R, 3, D> output = {range};

  //"test unipolar induced E from a dipole magnetic field"
  for (int j = range[th_].far_begin(); j < range[th_].far_end(); ++j) {
    for (int i = range[r_].far_begin(); i < range[r_].far_end(); ++i) {
      R cos_th = std::cos(grid[th_].absc(j, 0.5));
      E[r_]({i, j}) = std::exp(-2 * grid[r_].absc(i, 0)) * (1 - 3 * cos_th * cos_th);
      E[th_]({i, j}) = -std::exp(-2 * grid[r_].absc(i, 0.5)) *
                       std::sin(2 * grid[th_].absc(j, 0));
    }
  }

  SECTION("test dr  of \mathcal{E}_th") {
    output.reset();
    output[phi_] += diff_r(E,th_,1.0,b,e,false);

    R delta = grid[r_].delta();
    for (int j = 0; j < N+1; ++j) {
      for (int i = 0; i < N; ++i) {
        R x = 2 * std::sin(2*grid[th_].absc(j,0)) * std::exp(-2*grid[r_].absc(i,0));
        R corr = delta*delta;
        CAPTURE(i, j, corr, x, output[phi_]({i, j})-x);
        CHECK(output[phi_]({i, j}) == Approx(x).margin(corr));
      }
    }
  }

  SECTION("test d theta of \mathcal{E}_r ") {
    output.reset();
    output[phi_] += diff_th(E, r_, 1.0, b, e);

    R delta = grid[th_].delta();
    for (int j = 0; j < N+1; ++j) {
      for (int i = 0; i < N; ++i) {
        R x = 3 * std::sin(2 * grid[th_].absc(j, 0)) * std::exp(-2 * grid[r_].absc(i, 0));
        R corr = delta*delta;
        CAPTURE(i, j, corr, x, output[phi_]({i, j}) - x);
        CHECK(output[phi_]({i, j}) == Approx(x).margin(corr));
      }
    }
  }

}

TEST_CASE("test curlB", "[field]") {
  const int N = 128;
  const apt::Grid<R,D> grid = {{{0.0, std::log(30.0), N}, {0.0, std::acos(-1.0), N}}};
  constexpr int g = 1;

  auto range = apt::make_range<D>({}, {grid[0].dim(), grid[1].dim()}, g);

  std::optional<int> surf(0);
  std::optional<int> outer(N);
  std::optional<int> lo_axis(0);
  std::optional<int> hi_axis(N+1);

  const int b[D] = {0, 0};
  const int e[D] = {N, N+1};

  _helper<R,D> h(grid, 0, N, 0, N+1 );

  for (int i = 0; i < N; ++i)
    h.fr(i) = 1.0/std::exp(2.0*grid[r_].absc(i));

  Diff_in_CurlB<R, D> diffh_r{0, surf, outer, h, grid[r_].delta()};
  Diff_in_CurlB<R, D> diffh_th{1, surf, outer, h, grid[th_].delta()};

  Field<R, 3, D> B = {range};
  Field<R, 3, D> output = {range};

  //"test unipolar induced E from a dipole magnetic field"
  for (int j = range[th_].far_begin(); j < range[th_].far_end(); ++j) {
    for (int i = range[r_].far_begin(); i < range[r_].far_end(); ++i) {
      R cos_th = std::cos(grid[th_].absc(j, 0.0));
      B[r_]({i, j}) = std::exp(-2 * grid[r_].absc(i, 0.5)) *
                      (1 - 3 * cos_th * cos_th) *
                      std::sin(grid[th_].absc(j, 0.0));
      B[th_]({i, j}) = -std::exp(-2 * grid[r_].absc(i, 0.0)) *
                       std::sin(2 * grid[th_].absc(j, 0.5)) *
                       std::sin(grid[th_].absc(j, 0.5));
    }
  }

  SECTION("test dr  of \mathcal{B}_th") {
    output.reset();
    output[phi_] += diffh_r(B,th_,1.0,b,e);

    R delta = grid[r_].delta();
    for (int j = 0; j < N; ++j) {// NOTE j < N, not N+1, because cell N is beyond high axis
      for (int i = 0; i < N; ++i) {
        R x = 3 * std::sin(2*grid[th_].absc(j,0.5)) * std::exp(-4*grid[r_].absc(i,0.5)) * std::sin(grid[th_].absc(j,0.5));
        R corr = delta*delta;
        CAPTURE(i, j, corr, x, output[phi_]({i, j})-x);
        CHECK(output[phi_]({i, j}) == Approx(x).margin(corr));
      }
    }
  }

  SECTION("test d theta of \mathcal{B}_r ") {
    output.reset();
    output[phi_] += diffh_th(B, r_, 1.0, b, e);

    R delta = grid[th_].delta();
    for (int j = 0; j < N; ++j) { // NOTE upper bound is not N + 1
      for (int i = 0; i < N; ++i) {
        R x = 3 * std::sin(2 * grid[th_].absc(j, 0.5)) *
              std::exp(-4 * grid[r_].absc(i, 0.5)) *
              std::sin(grid[th_].absc(j, 0.5));
        R corr = delta*delta;
        CAPTURE(i, j, corr, x, output[phi_]({i, j}) - x);
        CHECK(output[phi_]({i, j}) == Approx(x).margin(corr));
      }
    }
  }

}
