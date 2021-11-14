#include <cassert>
#include <cmath>
#include <iostream>

#include "field/log_spherical_solver/updater.hpp"
#include "metric/log_spherical.hpp"
#include "pigeon.hpp"

// TODOL users may forget to sync value and state. Add another layer then
namespace particle {  // must have this namespace for now
template <typename T>
struct Specs {
  using value_type = T;
  static constexpr int Dim = 3;
  using state_type = long long;
  static_assert(8 * sizeof(state_type) >= 64);
};
}  // namespace particle

constexpr int DGrid = 2;
using real_t = float;
using real_j_t = float;
using real_export_t = float;
using ShapeF = particle::shapef_t<particle::shape::Cloud_In_Cell>;

using pgn = PIGEON<DGrid, real_t, particle::Specs, real_j_t, real_export_t>;

using Metric = metric::LogSpherical<real_t>;
using LogSphSolver_t = field::LogSphericalSolver<real_t, DGrid, real_j_t>;
constexpr real_t operator"" _r(long double x) noexcept {
  return static_cast<real_t>(x);
}

namespace {
// TODO already defined in apt/index.hpp, which is already included. But somehow
// these are not visible. Instead, vec_expression ops defineded in
// apt/numeric.hpp (in metric/log_spherical.hpp) is picked up and causing
// compile error.
template <int D>
constexpr apt::Index<D> operator+(apt::Index<D> ind,
                                  const apt::Longidx& l) noexcept {
  ind[l.dir()] += l.val();
  return ind;
}

constexpr apt::Longidx operator-(int a, apt::Longidx l) noexcept {
  l = a - l.val();
  return l;
}
}  // namespace

// TODO one issue with this is how to ensure that these values are set
struct global_variables {
  inline static real_t mu;
  inline static real_t Omega;
  inline static real_t spinup_time;

  inline static constexpr pgn::Grid_t supergrid = {
      {{0.0, std::log(30.0), 128}, {0.0, M_PI, 128}}};

  inline static constexpr int star_interior = 5;

  inline static constexpr int field_op_inv_precision = 4;
  inline static constexpr int guard = std::max(
      LogSphSolver_t::min_guard(field_op_inv_precision),
      (ShapeF::support() + 3) / 2);  // NOTE minimum number of guards of J
  // on one side is ( supp + 3 ) / 2
};
using gv = global_variables;

namespace {  // helper

real_t omega_spinup(real_t time) noexcept {
  return std::min<real_t>(time / gv::spinup_time, 1.0_r) * gv::Omega;
}

real_t B_r_star(real_t lnr, real_t theta, real_t, real_t time) noexcept {
  return gv::mu * 2.0_r * std::cos(theta) * std::exp(-3.0_r * lnr);
}
real_t B_theta_star(real_t lnr, real_t theta, real_t, real_t time) noexcept {
  return gv::mu * std::sin(theta) * std::exp(-3.0_r * lnr);
}
real_t B_phi_star(real_t lnr, real_t theta, real_t, real_t time) noexcept {
  return 0;
}

real_t E_r_star(real_t lnr, real_t theta, real_t, real_t time) noexcept {
  return gv::mu * omega_spinup(time) * std::exp(-2.0_r * lnr) *
         std::sin(theta) * std::sin(theta);
}
real_t E_theta_star(real_t lnr, real_t theta, real_t, real_t time) noexcept {
  return -gv::mu * omega_spinup(time) * std::exp(-2.0_r * lnr) *
         std::sin(2.0_r * theta);
}
real_t E_phi_star(real_t lnr, real_t theta, real_t, real_t time) noexcept {
  return 0;
}

template <typename T, typename F>
void axissymmetrize(
    field::Component<T, DGrid, false> comp,  // TODOL semantics on comp
    const F&
        f,  // has interface: void (*f)( real_t& val_guard, real_t& val_bulk ),
    pgn::Index_t Ib, pgn::Index_t Ie, bool is_upper) {
  static_assert(DGrid == 2);
  constexpr int AxisDir = 1;

  int mirror_sum = (comp.offset()[AxisDir] == MIDWAY) ? -1 : 0;
  mirror_sum += is_upper ? 2 * Ib[AxisDir] : 2 * (Ie[AxisDir] - 1);

  for (const auto& trI : apt::project_out(AxisDir, Ib, Ie)) {
    for (apt::Longidx n(AxisDir, Ib[AxisDir]); n < Ie[AxisDir]; ++n) {
      f(comp(trI + n), comp(trI + (mirror_sum - n)));
    }
  }
}

}  // namespace

namespace field {
class Solver : public pgn::FieldAction_t<Solver> {
 public:
  Solver(real_t four_pi, real_t alpha, int op_inv_precision, real_t surface,
         real_t outer)
      : m_solver(four_pi, alpha, op_inv_precision, surface, outer) {}

  void operator()(const Bundle_t& bundle) const override {
    m_solver(bundle.E, bundle.B, bundle.J, bundle.grid, bundle.cart,
             bundle.timestep, bundle.dt);
  }

 private:
  LogSphSolver_t m_solver;
};

// NOTE only implemented for UPPER boundary
struct RotatingConductor : public pgn::FieldAction_t<RotatingConductor> {
 private:
  apt::array<real_t (*)(real_t, real_t, real_t, real_t t), 3> m_E_cond;
  apt::array<real_t (*)(real_t, real_t, real_t, real_t t), 3> m_B_cond;

 public:
  auto& set_E_cond(
      apt::array<real_t (*)(real_t, real_t, real_t, real_t t), 3> E) {
    m_E_cond = E;
    return *this;
  }
  auto& set_B_cond(
      apt::array<real_t (*)(real_t, real_t, real_t, real_t t), 3> B) {
    m_B_cond = B;
    return *this;
  }

  virtual void operator()(const Bundle_t& bundle) const override {
    auto impose = [this, &grid = bundle.grid,
                   time = bundle.timestep * bundle.dt](
                      auto& F, const auto& F_cond, int comp) {
      const auto& ofs = F[comp].offset();
      const auto& m = F.mesh();
      const auto& ranges = this->ranges();

      for (int j = ranges[1].begin(); j < ranges[1].end(); ++j) {
        int li_j = m.linear_index(1, j);
        real_t th = grid[1].absc(j, ofs[1] * 0.5);
        for (int i = ranges[0].begin(); i < ranges[0].end(); ++i) {
          real_t lnr = grid[0].absc(i, ofs[0] * 0.5);
          int li = li_j + i * m.stride(0);
          F[comp][li] = F_cond[comp](lnr, th, 0, time);
        }
      }
    };
    for (int C = 0; C < 3; ++C) impose(bundle.E, m_E_cond, C);
    for (int C = 0; C < 3; ++C) impose(bundle.B, m_B_cond, C);
  }
};

// NOTE only implemented for UPPER boundary
struct DampingLayer : public pgn::FieldAction_t<DampingLayer> {
 private:
  apt::array<real_t (*)(real_t, real_t, real_t, real_t t), 3> m_E_bg;
  apt::array<real_t (*)(real_t, real_t, real_t, real_t t), 3> m_B_bg;
  real_t m_rate{};

  real_t m_r_damp_begin;  // only needed in m_profile
  real_t m_thickness;     // only_needed in m_profile
  real_t m_profile(real_t lnr) const {
    lnr = (std::exp(lnr) - m_r_damp_begin) / m_thickness;
    return 0.5 * lnr * lnr;
  }

 public:
  auto& set_E_background(
      apt::array<real_t (*)(real_t, real_t, real_t, real_t t), 3> E) {
    m_E_bg = E;
    return *this;
  }
  auto& set_B_background(
      apt::array<real_t (*)(real_t, real_t, real_t, real_t t), 3> B) {
    m_B_bg = B;
    return *this;
  }
  auto& set_damping_rate(real_t rate) {
    m_rate = rate;
    return *this;
  }

  auto& set_damping_begin(real_t r_damp_b) {
    m_r_damp_begin = r_damp_b;
    return *this;
  }

  auto& set_damping_thickness(real_t thickness) {
    m_thickness = thickness;
    return *this;
  }

  void operator()(const Bundle_t& bundle) const override {
    auto impose = [this, &grid = bundle.grid, dt = bundle.dt](
                      auto& F, const auto& F_bg, int comp) {
      const auto& ofs = F[comp].offset();
      const auto& m = F.mesh();
      const auto& ranges = this->ranges();

      for (int j = ranges[1].begin(); j < ranges[1].end(); ++j) {
        int li_j = m.linear_index(1, j);
        real_t th = grid[1].absc(j, ofs[1] * 0.5);
        for (int i = ranges[0].begin(); i < ranges[0].end(); ++i) {
          int li = li_j + i * m.stride(0);
          real_t lnr = grid[0].absc(i, ofs[0] * 0.5);
          real_t lambda = 1.0 - m_rate * dt * m_profile(lnr);
          real_t f_bg =
              F_bg[comp](lnr, th, 0.0,
                         0.0);  // time = 0.0, so damp to the initial condition
          F[comp][li] = (F[comp][li] - f_bg) * lambda + f_bg;
        }
      }
    };
    for (int C = 0; C < 3; ++C) impose(bundle.E, m_E_bg, C);
    for (int C = 0; C < 3; ++C) impose(bundle.B, m_B_bg, C);
  }
};

struct Axissymmetric : public pgn::FieldAction_t<Axissymmetric> {
 private:
  bool m_is_upper_axis = false;

 public:
  auto& is_upper_axis(bool x) {
    m_is_upper_axis = x;
    return *this;
  }

  void operator()(const Bundle_t& bundle) const override {
    // NOTE Guard cells values are needed when interpolating E and B
    // E_theta, B_r, B_phi are on the axis. All but B_r should be set to zero
    auto assign = [](real_t& v_g, real_t& v_b) noexcept { v_g = v_b; };
    auto neg_assign = [](real_t& v_g, real_t& v_b) noexcept {
      v_g = (&v_g == &v_b) ? 0.0_r : -v_b;
    };
    const auto& ranges = this->ranges();
    const auto begins = apt::range::begin(ranges);
    const auto ends = apt::range::end(ranges);

    // MIDWAY in AxisDir
    axissymmetrize(bundle.E[0], assign, begins, ends, m_is_upper_axis);
    axissymmetrize(bundle.E[2], neg_assign, begins, ends, m_is_upper_axis);
    axissymmetrize(bundle.B[1], neg_assign, begins, ends, m_is_upper_axis);

    // INSITU in AxisDir
    axissymmetrize(bundle.E[1], neg_assign, begins, ends, m_is_upper_axis);
    axissymmetrize(bundle.B[0], assign, begins, ends, m_is_upper_axis);
    axissymmetrize(bundle.B[2], neg_assign, begins, ends, m_is_upper_axis);
  }
};
}  // namespace field

auto set_up_particle_properties() {
  particle::map<particle::Properties> properties;
  return properties;
}

int main(int argc, char** argv) {
  const auto args = pic::parse_args(argc, argv);
  assert(args.config_file);

  const apt::array<int, DGrid> dims = {1, 1};
  const apt::array<bool, DGrid> periodic = {false, false};
  gv::Omega = 1.0 / 6.0;

  auto conf = pgn::ConfFile_t::load(*args.config_file);

  const std::string datadir_prefix =
      conf["datadir_prefix"].as_or<std::string>("../Data/");
  const std::string project_name =
      conf["project_name"].as_or<std::string>("Unnamed");
  const int n_timesteps = conf["total_timesteps"].as_or<int>(100);

  const real_t dt = conf["dt"].as<real_t>();
  gv::mu = conf["gamma0"].as<real_t>() * std::pow(gv::Omega, -2.0);
  const real_t wpic2 = 2 * gv::Omega * gv::mu / conf["Np"].as<real_t>();

  const real_t gamma_fd = conf["pairs"]["gamma_fd"].as<real_t>();
  const real_t E_ph = conf["pairs"]["E_ph"].as<real_t>();
  const real_t gravity_strength = conf["forces"]["gravity"].as<real_t>();
  const real_t landau0_B_thr =
      conf["forces"]["landau0_ratio"].as<real_t>() * gv::mu;
  const std::array<real_t, 2> mfp = {
      conf["pairs"]["photon"]["mfp"][0].as<real_t>(),
      conf["pairs"]["photon"]["mfp"][1]
          .as<real_t>()};  // TODO this should use toml array parsing

  const real_t v_th = conf["atmosphere"]["v_th"].as<real_t>();
  const real_t atm_x = conf["atmosphere"]["multiplicity"].as<real_t>();

  const real_t damping_thickness = conf["damping"]["thickness"].as<real_t>();
  const real_t damping_rate = conf["damping"]["rate"].as<real_t>();
  gv::spinup_time = conf["spinup_time"].as<real_t>();

  const real_t r_e = wpic2 * apt::dV(gv::supergrid) / (4.0 * M_PI);

  pgn::SimulationBuilder_t builder(args);

  builder.initialize_this_run_dir(datadir_prefix, project_name)
      .create_cartesian_topology(dims, periodic);

  // TODO
  auto properties = set_up_particle_properties();
  builder.set_particle_properties(properties);

  {  // set up field actions
    builder
        .add_field_action<field::Solver>(
            4.0 * M_PI * r_e, 1.0, gv::field_op_inv_precision,
            gv::supergrid[0].absc(gv::star_interior), gv::supergrid[0].upper())
        .set_name("LogSphericalSolver")
        .set_range(0, {gv::star_interior, gv::supergrid[0].dim(), gv::guard})
        .set_range(1, {0, gv::supergrid[1].dim() + 1, gv::guard});

    builder.add_field_action<field::RotatingConductor>()
        .set_E_cond({E_r_star, E_theta_star, E_phi_star})
        .set_B_cond({B_r_star, B_theta_star, B_phi_star})
        .set_name("RotatingConductor")
        .set_range(0, {0, gv::star_interior})
        .set_range(1, {0, gv::supergrid[1].dim() + 1});

    const int damping_begin_index = gv::supergrid[0].csba(
        std::log(std::exp(gv::supergrid[0].upper()) - damping_thickness));
    builder.add_field_action<field::DampingLayer>()
        .set_E_background({E_r_star, E_theta_star, E_phi_star})
        .set_B_background({B_r_star, B_theta_star, B_phi_star})
        .set_damping_rate(damping_rate)
        .set_damping_begin(std::exp(gv::supergrid[0].upper()) -
                           damping_thickness)
        .set_damping_thickness(damping_thickness)
        .set_name("DampingLayer")
        .set_range(0, {damping_begin_index, gv::supergrid[0].dim() + gv::guard})
        .set_range(1, {0, gv::supergrid[1].dim() + 1});

    builder.add_field_action<field::Axissymmetric>()
        .is_upper_axis(false)
        .set_name("AxissymmetrizeEBLower")
        .set_range(0, {0, gv::supergrid[0].dim()})
        .set_range(
            1, {-gv::guard, 1});  // NOTE 1 so as to set values right on axis

    builder.add_field_action<field::Axissymmetric>()
        .is_upper_axis(true)
        .set_name("AxissymmetrizeEBHigher")
        .set_range(0, {0, gv::supergrid[0].dim()})
        .set_range(
            1, {gv::supergrid[1].dim(), gv::supergrid[1].dim() + gv::guard});
  }

  auto& sim = builder.build();

  // TODO having to have user call this is a bit error_prone
  const auto init_ts = sim.initial_timestep();
  if (mpi::world.rank() == 0) std::cout << "Launch" << std::endl;
  for (int ts = init_ts; ts < init_ts + n_timesteps; ++ts) {
    sim.evolve(ts, dt);
  }

  return 0;
}
