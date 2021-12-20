#include <cassert>
#include <cmath>
#include <iostream>

#include "field/log_spherical_solver/updater.hpp"
#include "metric/log_spherical.hpp"
#include "msh/mesh_shape_interplay.hpp"
#include "particle/annihilation.hpp"
#include "particle/forces.hpp"
#include "particle/scattering.hpp"
#include "particle/updater.hpp"
#include "pic/forces/gravity.hpp"
#include "pic/forces/landau0.hpp"
#include "pigeon.hpp"

// TODOL users may forget to sync value and state. Add another layer then
namespace particle {  // must have this namespace for now
template <typename T>
struct Specs {
  using value_type = T;
  static constexpr int Dim = 3;
  using state_type =
      apt::copy_cvref_t<T, long long>;  // NOTE needed to compile e.g. the ptc
                                        // loop in f_count() in Atmosphere,
                                        // where const long long is needed.
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
using ParticleUpdater_t =
    particle::Updater<DGrid, real_t, particle::Specs, ShapeF, real_j_t>;

// define convenient literal operators
constexpr real_t operator"" _r(long double x) noexcept {
  return static_cast<real_t>(x);
}

constexpr real_t operator"" _deg(long double x) noexcept {
  return static_cast<real_t>(x * M_PI / 180);
}

// TODO one issue with this is how to ensure that these values are set
struct global_variables {
  inline static real_t mu;
  inline static real_t Omega;
  inline static real_t wpic2;
  inline static real_t spinup_time;

  inline static constexpr pgn::Grid_t supergrid = {
      {{0.0, std::log(30.0), 128}, {0.0, M_PI, 128}}};

  inline static constexpr int star_interior = 5;

  inline static constexpr int field_op_inv_precision = 4;
  inline static constexpr int guard = std::max(
      LogSphSolver_t::min_guard(field_op_inv_precision),
      (ShapeF::support() + 3) / 2);  // NOTE minimum number of guards of J
  // on one side is ( supp + 3 ) / 2

  /** The following are some custom data
   */
  inline static particle::map<real_t> N_scat{};
  inline static pgn::Field<4> pc_counter{};
  inline static pgn::Field<4> em_counter{};
  inline static real_t cumulative_time{};

  inline static std::optional<particle::map<pgn::JField>>
      Jsp;  // current by species
  inline static std::optional<pgn::Field<1>> skin_depth;
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
void axisymmetrize(
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

struct Axisymmetric : public pgn::FieldAction_t<Axisymmetric> {
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
    axisymmetrize(bundle.E[0], assign, begins, ends, m_is_upper_axis);
    axisymmetrize(bundle.E[2], neg_assign, begins, ends, m_is_upper_axis);
    axisymmetrize(bundle.B[1], neg_assign, begins, ends, m_is_upper_axis);

    // INSITU in AxisDir
    axisymmetrize(bundle.E[1], neg_assign, begins, ends, m_is_upper_axis);
    axisymmetrize(bundle.B[0], assign, begins, ends, m_is_upper_axis);
    axisymmetrize(bundle.B[2], neg_assign, begins, ends, m_is_upper_axis);
  }
};
}  // namespace field

auto set_up_particle_properties(real_t gravity_strength, real_t landau0_B_thr,
                                real_t gamma_fd, real_t gamma_off,
                                real_t Ndot_fd, real_t E_ph,
                                real_t magnetic_convert_B_thr,
                                const std::array<real_t, 2>& mfp) {
  using namespace particle;
  map<Properties> properties;
  {
    properties.insert(species::electron, {1.0, -1.0, "electron", "el"});
    properties.insert(species::positron, {1.0, 1.0, "positron", "po"});
    properties.insert(species::ion, {5.0, 1.0, "ion", "io"});
    properties.insert(species::photon, {0, 0, "photon", "ph"});
  }

  {
    using Force = Force<real_t, Specs>;
    constexpr auto* lorentz = force::template lorentz<real_t, Specs, vParticle>;
    constexpr auto* landau0 = force::landau0<real_t, Specs, vParticle>;
    constexpr auto* gravity = force::gravity<real_t, Specs, vParticle>;

    if (properties.has(species::electron)) {
      auto sp = species::electron;
      Force force;
      const auto& prop = properties[sp];

      force.add(lorentz, prop.charge_x / prop.mass_x);
      force.add(gravity, gravity_strength);
      force.add(landau0, landau0_B_thr);

      force.Register(sp);
    }
    if (properties.has(species::positron)) {
      auto sp = species::positron;
      Force force;
      const auto& prop = properties[sp];

      force.add(lorentz, prop.charge_x / prop.mass_x);
      force.add(gravity, gravity_strength);
      force.add(landau0, landau0_B_thr);

      force.Register(sp);
    }
    if (properties.has(species::ion)) {
      auto sp = species::ion;
      Force force;
      const auto& prop = properties[sp];

      force.add(lorentz, prop.charge_x / prop.mass_x);
      force.add(gravity, gravity_strength);
      // force.add( landau0, landau0_B_thr );

      force.Register(sp);
    }
  }

  using Ptc_t = typename pgn::ParticleArray_t::particle_type;

  static auto flagger = [](flagbits parent_bits, species) noexcept {
    parent_bits[flag::secondary] = true;
    return parent_bits;
  };

  {
    Scat<real_t, Specs> ep_scat;

    ep_scat.eligs.push_back(
        [](const Ptc_t& ptc) { return ptc.q(0) < std::log(9.0_r); });

    scat::CurvatureRadiation<real_t, Specs>::gamma_fd = gamma_fd;
    scat::CurvatureRadiation<real_t, Specs>::gamma_off = gamma_off;
    scat::CurvatureRadiation<real_t, Specs>::Ndot_fd = Ndot_fd;
    scat::CurvatureRadiation<real_t, Specs>::E_ph = E_ph;
    ep_scat.channels.push_back(scat::CurvatureRadiation<real_t, Specs>::test);

    if (properties.has(species::photon))
      ep_scat.impl = [](auto itr, auto& p, real_t t) {
        return scat::RadiationFromCharges<false, real_t, Specs>(itr, p, t,
                                                                flagger);
      };
    else
      ep_scat.impl = [](auto itr, auto& p, real_t t) {
        return scat::RadiationFromCharges<true, real_t, Specs>(itr, p, t,
                                                               flagger);
      };

    if (properties.has(species::electron) &&
        properties.has(species::positron)) {
      ep_scat.Register(species::electron);
      ep_scat.Register(species::positron);
    }
  }

  if (properties.has(species::photon)) {
    ::particle::Scat<real_t, Specs> photon_scat;
    photon_scat.eligs.push_back([](const Ptc_t& ptc) { return true; });
    scat::MagneticConvert<real_t, Specs>::B_thr = magnetic_convert_B_thr;
    scat::MagneticConvert<real_t, Specs>::mfp = mfp[0];
    photon_scat.channels.push_back(scat::MagneticConvert<real_t, Specs>::test);

    scat::TwoPhotonCollide<real_t, Specs>::mfp = mfp[1];
    photon_scat.channels.push_back(scat::TwoPhotonCollide<real_t, Specs>::test);

    photon_scat.impl = [](auto itr, auto& p, real_t t) {
      return scat::PhotonPairProduction<real_t, Specs>(itr, p, t, flagger);
    };

    photon_scat.Register(species::photon);
  }
  return properties;
}

namespace particle {
struct MainUpdater : public pgn::ParticleAction_t<MainUpdater> {
 public:
  auto& set_update_q(ParticleUpdater_t::UpdateQ_t update_q) {
    m_updater.set_update_q(update_q);
    return *this;
  }

  auto& set_ignore_current(std::unordered_set<species> sps) {
    m_updater.set_ignore_current(std::move(sps));
    return *this;
  }

  auto& set_export_schedule(pgn::ExportSchedule sch) {
    m_sch = std::move(sch);
    return *this;
  }

  void operator()(const Bundle_t& bd) const override {
    auto run_updater = [&](auto& ptcs, auto& J) {
      m_updater(ptcs, J, &bd.new_ptc_buf, bd.properties, bd.E, bd.B, bd.grid,
                bd.dt, bd.timestep, bd.rng);
    };

    if (gv::Jsp and m_sch.is_do(bd.timestep)) {
      auto& Jsp = *gv::Jsp;
      // store J by species separately for data export
      for (auto sp : bd.particles) {
        map<pgn::ParticleArray_t> ptcs_sp;
        ptcs_sp.insert(sp);
        std::swap(ptcs_sp[sp], bd.particles[sp]);
        Jsp[sp] = pgn::JField(bd.J.mesh());
        run_updater(ptcs_sp, Jsp[sp]);
        std::swap(ptcs_sp[sp], bd.particles[sp]);

        for (int C = 0; C < 3; ++C) {
          for (int i = 0; i < bd.J.mesh().linear_size(); ++i)
            bd.J[C][i] += Jsp[sp][C][i];
        }
      }
    } else {
      run_updater(bd.particles, bd.J);
    }
  }

 private:
  ParticleUpdater_t m_updater;
  pgn::ExportSchedule m_sch;
};

struct Atmosphere : public pgn::ParticleAction_t<Atmosphere> {
 private:
  mutable pgn::Field<1> m_count_n;
  mutable pgn::Field<1> m_count_p;
  int m_n = 0;  // normal direction
  species m_posion = species::ion;
  species m_negaon = species::electron;
  real_t m_v_th = 0.0;
  real_t m_N_atm = 0.0;
  real_t m_min_frac = 1e-6;  // over fracs larger than this will be injected
  real_t (*m_omega_t)(real_t time) = nullptr;

 public:
  auto& set_thermal_velocity(real_t v) {
    m_v_th = v;
    return *this;
  }
  auto& set_number_in_atmosphere(real_t N) {
    m_N_atm = N;
    return *this;
  }
  auto& set_minimal_fraction(real_t x) {
    m_min_frac = x;
    return *this;
  }
  auto& set_normal_direction(int n) {
    m_n = n;
    return *this;
  }
  auto& set_omega_t(real_t (*omega_t)(real_t)) {
    m_omega_t = omega_t;
    return *this;
  }
  auto& set_positive_charge(species sp) {
    m_posion = sp;
    return *this;
  }
  auto& set_negative_charge(species sp) {
    m_negaon = sp;
    return *this;
  }

  void operator()(const Bundle_t& bundle) const override {
    const auto& grid = bundle.grid;
    const auto& ens = bundle.ens;
    const auto& B = bundle.B;
    auto& particles = bundle.particles;
    auto& rng = bundle.rng;
    auto timestep = bundle.timestep;
    auto dt = bundle.dt;

    const auto& ranges = this->ranges();
    const auto begins = apt::range::begin(ranges);
    const auto ends = apt::range::end(ranges);

    if (ranges[m_n].end() <= ranges[m_n].begin()) return;

    m_count_n.resize({apt::make_range(begins, ends, 0)});
    m_count_p.resize({apt::make_range(begins, ends, 0)});

    apt::array<real_t, DGrid> lb;
    apt::array<real_t, DGrid> ub;
    for (int i = 0; i < DGrid; ++i) {
      lb[i] = grid[i].absc(ranges[i].begin(), 0.0);
      ub[i] = grid[i].absc(ranges[i].end(), 0.0);
    }
    auto is_in = [&lb, &ub](const auto& q) {
      for (int i = 0; i < DGrid; ++i) {
        if (q[i] < lb[i] || q[i] >= ub[i]) return false;
      }
      return true;
    };

    auto f_count = [&lb, &ub, &grid, is_in](auto& count, const auto& ptcs) {
      count.reset();
      for (const auto& x : ptcs) {
        if (!x.is(flag::exist) || !is_in(x.q())) continue;
        pgn::Index_t idx;
        for (int i = 0; i < DGrid; ++i) {
          // NOTE used grid.lower instead of lb, important
          idx[i] = grid[i].csba(x.q(i));
        }
        count[0](idx) += x.frac();  // add by fraction
      }
    };

    f_count(m_count_n, particles[m_negaon]);
    f_count(m_count_p, particles[m_posion]);

    {  // parallelize TODO optimize
      int rank_inj = timestep % ens.size();
      ens.intra.template reduce<true>(mpi::by::SUM, rank_inj,
                                      m_count_n[0].data().data(),
                                      m_count_n[0].data().size());
      ens.intra.template reduce<true>(mpi::by::SUM, rank_inj,
                                      m_count_p[0].data().data(),
                                      m_count_p[0].data().size());
      if (ens.intra.rank() != rank_inj) return;
    }

    auto itr_po = std::back_inserter(particles[m_posion]);
    auto itr_ne = std::back_inserter(particles[m_negaon]);

    for (const auto& I : apt::Block(begins, ends)) {
      auto N_pairs = std::min(m_count_n[0](I), m_count_p[0](I));
      pgn::Vec3 q{};
      for (int i = 0; i < DGrid; ++i) q[i] = grid[i].absc(I[i], 0.5);

      pgn::Vec3 nB{};
      {  // make nB centered in the cell
        const auto& m = B.mesh();
        auto li = m.linear_index(I);
        if constexpr (DGrid == 2) {
          nB[0] = 0.5_r * (B[0][li] + B[0][li + m.stride(1)]);
          nB[1] = 0.5_r * (B[1][li] + B[1][li + m.stride(0)]);
          nB[2] = 0.25_r *
                  (B[2][li] + B[2][li + m.stride(0)] + B[2][li + m.stride(1)] +
                   B[2][li + m.stride(0) + m.stride(1)]);
        } else if (DGrid == 3) {
          nB[0] = 0.25_r *
                  (B[0][li] + B[0][li + m.stride(1)] + B[0][li + m.stride(2)] +
                   B[0][li + m.stride(1) + m.stride(2)]);
          nB[1] = 0.25_r *
                  (B[1][li] + B[1][li + m.stride(2)] + B[1][li + m.stride(0)] +
                   B[1][li + m.stride(2) + m.stride(0)]);
          nB[2] = 0.25_r *
                  (B[2][li] + B[2][li + m.stride(0)] + B[2][li + m.stride(1)] +
                   B[2][li + m.stride(0) + m.stride(1)]);
        }
        if (apt::abs(nB) == 0.0_r)
          nB = {1.0_r, 0.0_r, 0.0_r};  // use radial direction as default
        else
          nB /= apt::abs(nB);
      }

      pgn::Vec3 p{};
      p[2] = m_omega_t(timestep * dt) * std::exp(q[0]) *
             std::sin(q[1]);  // corotating

      // replenish
      real_t quota = m_N_atm * std::sin(q[1]) - N_pairs;
      while (quota > m_min_frac) {
        auto q_ptc = q;
        real_t frac = std::min(1.0_r, quota);
        quota -= 1.0_r;

        for (int i = 0; i < DGrid; ++i) {
          if (m_n == i)
            q_ptc[i] += grid[i].delta() * rng.uniform(-0.5, 0.0);
          else
            q_ptc[i] += grid[i].delta() * rng.uniform(-0.5, 0.5);
        }
        auto p_ptc = p;
        p_ptc += nB * rng.gaussian(0.0, m_v_th);
        *(itr_ne++) = pgn::Particle(q_ptc, p_ptc, frac, m_negaon);
        *(itr_po++) =
            pgn::Particle(std::move(q_ptc), std::move(p_ptc), frac, m_posion);
      }
    }
  }
};

struct Axisymmetric : public pgn::ParticleAction_t<Axisymmetric> {
 private:
  bool _is_upper_axis = false;

 public:
  auto& is_upper_axis(bool x) {
    _is_upper_axis = x;
    return *this;
  }

  void operator()(const Bundle_t& bundle) const override {
    const auto& grid = bundle.grid;
    auto& J = bundle.J;
    auto add_assign = [](real_j_t& a, real_j_t& b) noexcept {
      a += b;
      b = a;
    };

    auto sub_assign = [](real_j_t& a, real_j_t& b) noexcept {
      a -= b;
      b = -a;
    };

    const auto& ranges = this->ranges();
    const auto begins = apt::range::begin(ranges);
    const auto ends = apt::range::end(ranges);
    // MIDWAY in AxisDir
    axisymmetrize(J[0], add_assign, begins, ends, _is_upper_axis);
    axisymmetrize(J[2], sub_assign, begins, ends, _is_upper_axis);
    // INSITU in AxisDir
    axisymmetrize(J[1], sub_assign, begins, ends, _is_upper_axis);
  }
};

struct NewPtcAnalyzer : public pgn::ParticleAction_t<NewPtcAnalyzer> {
  void operator()(const Bundle_t& bd) const override {
    // Put particles where they belong after scattering
    for (auto& ptc : bd.new_ptc_buf) {
      if (not ptc.is(flag::exist)) continue;
      const auto this_sp = ptc.template get<species>();

      // FIXME this relys on all particles in this buffer that carry
      // flag::secondary are new particles. One possible exception is
      // migration. So this action must be done before migration.
      if (ptc.is(flag::secondary)) {
        // log scattering events
        gv::N_scat[this_sp] += ptc.frac();

        // log pair creation events
        if (species::electron == this_sp) {
          pgn::Index_t I;  // domain index, not the global index
          for (int j = 0; j < DGrid; ++j) I[j] = bd.grid[j].csba(ptc.q(j));
          gv::pc_counter[0](I) += ptc.frac();
          for (int j = 0; j < 3; ++j)
            gv::pc_counter[j + 1](I) +=
                2.0 * ptc.frac() * ptc.p(j);  // 2 because of positron
        } else if (species::photon == this_sp) {
          pgn::Index_t I;  // domain index, not the global index
          for (int j = 0; j < DGrid; ++j) I[j] = bd.grid[j].csba(ptc.q(j));
          gv::em_counter[0](I) += ptc.frac();
          for (int j = 0; j < 3; ++j)
            gv::em_counter[j + 1](I) += ptc.frac() * ptc.p(j);
        }
      }

      bd.particles[this_sp].push_back(std::move(ptc));
    }
    bd.new_ptc_buf.resize(0);
    gv::cumulative_time += bd.dt;
  }
};

struct Escaping : public pgn::ParticleAction_t<Escaping> {
 private:
  int m_n = 0;  // normal direction

 public:
  void operator()(const Bundle_t& bundle) const override {
    const auto& ranges = this->ranges();
    if (apt::range::is_empty(ranges)) return;

    auto& particles = bundle.particles;
    const auto& grid = bundle.grid;

    for (auto sp : particles) {
      for (auto ptc : particles[sp]) {  // TODOL semantics
        // NOTE: this implementation assumes that there is never particles
        // traveling backwards. If not, the flag::ignore_force must be unset
        // when such particles leave the layer, and one must be careful that
        // the ignore_force delimiter must not coincide with the lower bound
        // of a patch.
        if (ptc.q(m_n) >= grid[m_n].absc(ranges[m_n].begin()) and
            ptc.q(m_n) < grid[m_n].absc(ranges[m_n].end())) {
          if (not ptc.is(flag::ignore_force)) {
            // turn momentum to radial
            auto p = apt::abs(ptc.p());
            ptc.p(1) = 0;
            ptc.p(2) = 0;
            ptc.p(m_n) = p;
            ptc.set(flag::ignore_force);
          }
        }
      }
    }
  }
};

class Annihilator : public pgn::ParticleAction_t<Annihilator> {
 public:
  void operator()(const Bundle_t& bundle) const override {
    const auto& ranges = this->ranges();
    if (apt::range::is_empty(ranges)) return;

    const auto& grid = bundle.grid;
    const auto& ens = bundle.ens;
    auto& particles = bundle.particles;
    const auto& properties = bundle.properties;
    auto& J = bundle.J;
    auto dt = bundle.dt;
    auto timestep = bundle.timestep;

    if (m_sch.is_do(timestep)) {
      apt::array<apt::array<real_t, 2>, DGrid> bounds;
      for (int i = 0; i < DGrid; ++i) {
        // NOTE guard cells or margins of range not included.
        bounds[i][0] = grid[i].absc(ranges[i].begin());
        bounds[i][1] = grid[i].absc(ranges[i].end());
      }

      particle::annihilate(particles[species::electron],
                           particles[species::positron], J,
                           properties[species::electron].charge_x,
                           properties[species::positron].charge_x, grid,
                           ens.intra, dt, ShapeF(), m_policy, bounds);
    }
  }

  auto& schedule() { return m_sch; }

 private:
  pgn::Schedule m_sch;

  // policy returns number of pairs to be annihilated
  static real_t m_policy(real_t num_e, real_t num_p) noexcept {
    return std::min(num_e, num_p) / 2;
  }
};

class NoPhotonZone : public pgn::ParticleAction_t<NoPhotonZone> {
 private:
  int m_interval = 100;

 public:
  auto& set_interval(int interval) noexcept {
    m_interval = interval;
    return *this;
  }

  void operator()(const Bundle_t& bundle) const override {
    const auto& ranges = this->ranges();
    auto timestep = bundle.timestep;
    const auto& grid = bundle.grid;
    auto& particles = bundle.particles;

    if ((timestep % m_interval != 0) or apt::range::is_empty(ranges)) return;

    apt::array<apt::array<real_t, 2>, DGrid> bounds;
    for (int i = 0; i < DGrid; ++i) {
      bounds[i][0] = grid[i].absc(ranges[i].begin());
      bounds[i][1] = grid[i].absc(ranges[i].end());
    }

    auto in_zone = [&bounds](const auto& q) noexcept {
      for (int i = 0; i < DGrid; ++i) {
        if (q[i] < bounds[i][0] or q[i] >= bounds[i][1]) return false;
      }
      return true;
    };

    for (auto ptc : particles[species::photon]) {  // TODOL semantics
      if (ptc.is(flag::exist) and in_zone(ptc.q())) {
        ptc.reset(flag::exist);
      }
    }
  }
};

}  // namespace particle

struct InitialCondition
    : public pgn::InitialConditionAction_t<InitialCondition> {
  void operator()(const Bundle_t& bundle) const override {
    const auto& ranges = this->ranges();
    const auto& grid = bundle.grid;
    auto& B = bundle.B;
    for (const auto& I :
         apt::Block(apt::range::begin(ranges), apt::range::end(ranges))) {
      B[0](I) = B_r_star(grid[0].absc(I[0], 0.5 * B[0].offset()[0]),
                         grid[1].absc(I[1], 0.5 * B[0].offset()[1]), 0, 0);
      B[1](I) = B_theta_star(grid[0].absc(I[0], 0.5 * B[1].offset()[0]),
                             grid[1].absc(I[1], 0.5 * B[1].offset()[1]), 0, 0);
    }
  }
};

struct PostResume : public pgn::PostResumeAction_t<PostResume> {
  void operator()(const Bundle_t& bundle) const override {
    // not implemented
  }
};

namespace io {
void do_prior_export(const pgn::ExportBundle_t& bd) {
  {  // pair creation counter
    auto& pc = gv::pc_counter;
    for (int i = 0; i < 4; ++i)
      bd.ens.reduce_to_chief(mpi::by::SUM, pc[i].data().data(),
                             pc[i].data().size());
  }
  {  // photon emission counter
    auto& em = gv::em_counter;
    for (int i = 0; i < 4; ++i)
      bd.ens.reduce_to_chief(mpi::by::SUM, em[i].data().data(),
                             em[i].data().size());
  }
  if (gv::Jsp) {
    auto& Jsp = *gv::Jsp;
    for (auto s : Jsp) {
      for (int i = 0; i < 3; ++i)
        bd.ens.reduce_to_chief(mpi::by::SUM, Jsp[s][i].data().data(),
                               Jsp[s][i].data().size());
    }
    if (bd.cart_opt) {
      for (auto s : Jsp) {
        field::merge_sync_guard_cells(Jsp[s], *bd.cart_opt);
      }
    }
  }

  if (gv::skin_depth) {  // skin depth
    auto& skd = *gv::skin_depth;
    skd = {bd.J.mesh()};
    skd.reset();

    for (auto sp : bd.particles) {
      auto q2m = bd.properties[sp].charge_x * bd.properties[sp].charge_x /
                 bd.properties[sp].mass_x;
      for (const auto& ptc : bd.particles[sp]) {
        if (!ptc.is(particle::flag::exist)) continue;
        pgn::Index_t I;
        for (int i = 0; i < DGrid; ++i) I[i] = bd.grid[i].csba(ptc.q(i));
        skd[0](I) += q2m * ptc.frac();
      }
    }
    bd.ens.reduce_to_chief(mpi::by::SUM, skd[0].data().data(),
                           skd[0].data().size());
    if (bd.ens.is_chief()) {
      for (const auto& I : apt::Block(apt::range::begin(skd.mesh().range()),
                                      apt::range::end(skd.mesh().range()))) {
        real_t r = bd.grid[0].absc(I[0], 0.5);
        real_t theta = bd.grid[1].absc(I[1], 0.5);
        real_t h = Metric::h<2>(r, theta) /
                   (gv::wpic2 * bd.grid[0].delta() * bd.grid[0].delta());
        auto& v = skd[0](I);
        v = std::sqrt(h / v);
      }
    }
  }
}

void do_post_export(const pgn::ExportBundle_t&) {
  auto& pc = gv::pc_counter;
  for (int i = 0; i < 4; ++i)
    std::fill(pc[i].data().begin(), pc[i].data().end(), 0);
  auto& em = gv::em_counter;
  for (int i = 0; i < 4; ++i)
    std::fill(em[i].data().begin(), em[i].data().end(), 0);
  gv::cumulative_time = 0;
  if (gv::Jsp) {
    // clear Jsp to save some space
    for (auto sp : *gv::Jsp) (*gv::Jsp)[sp] = {};
  }
  if (gv::skin_depth) *gv::skin_depth = {};
}

constexpr auto I2std(const pgn::Index_t& I) {
  apt::array<real_t, DGrid> res;
  for (int i = 0; i < DGrid; ++i) res[i] = I[i] + 0.5;  // interpolate to MIDWAY
  return res;
}

template <int F>
apt::array<real_t, 3> field_self(pgn::Index_t I, const pgn::Grid_t& grid,
                                 const pgn::Field<3>& E, const pgn::Field<3>& B,
                                 const pgn::JField& J) {
  if constexpr (F == 0) {
    return msh::interpolate(E, I2std(I), ShapeF());
  } else if (F == 1) {
    return msh::interpolate(B, I2std(I), ShapeF());
  } else if (F == 2) {
    return msh::interpolate(J, I2std(I), ShapeF());
  } else {
    static_assert(F < 3);
  }
}

apt::array<real_t, 3> Br_alone(pgn::Index_t I, const pgn::Grid_t& grid,
                               const pgn::Field<3>& E, const pgn::Field<3>& B,
                               const pgn::JField& J) {
  return {msh::interpolate(B[0], I2std(I), ShapeF()), 0, 0};
}

constexpr int POW(int B, int E) {
  if (E == 0)
    return 1;
  else
    return B * POW(B, E - 1);
}

void average_when_downsampled(pgn::IOField& fds, const pgn::IOGrid& grid,
                              int num_comps, const mpi::CartComm&) {
  int downsample_ratio =
      static_cast<int>(grid[0].delta() / gv::supergrid[0].delta() + 0.5);
  int factor = POW(downsample_ratio, DGrid);
  for (int i = 0; i < num_comps; ++i) {
    for (auto& x : fds[i].data()) x /= factor;
  }
}

void average_and_divide_flux_by_area(pgn::IOField& fds, const pgn::IOGrid& grid,
                                     int num_comps, const mpi::CartComm& cart) {
  average_when_downsampled(fds, grid, num_comps, cart);

  using Metric = metric::LogSpherical<real_export_t>;
  // define a function pointer.
  real_export_t (*hh_func)(real_export_t, real_export_t, real_export_t) =
      nullptr;
  apt::array<real_export_t, 3> q{};

  for (int comp = 0; comp < num_comps; ++comp) {
    const auto& ofs = fds[comp].offset();
    switch (comp) {
      case 0:
        hh_func = Metric::template hh<0>;
        break;
      case 1:
        hh_func = Metric::template hh<1>;
        break;
      case 2:
        hh_func = Metric::template hh<2>;
        break;
    }

    for (const auto& I : apt::Block(apt::range::begin(fds.mesh().range()),
                                    apt::range::end(fds.mesh().range()))) {
      for (int i = 0; i < DGrid; ++i) q[i] = grid[i].absc(I[i], 0.5 * ofs[i]);
      auto hh = hh_func(q[0], q[1], q[2]);
      if (std::abs(hh) > 1e-12)
        fds[comp](I) /= hh;
      else
        fds[comp](I) = 0.0;
    }
  }
}

apt::array<real_t, 3> pair_creation_rate(pgn::Index_t I,
                                         const pgn::Grid_t& grid,
                                         const pgn::Field<3>&,
                                         const pgn::Field<3>&,
                                         const pgn::JField&) {
  auto x = gv::pc_counter[0](I);
  return {x / gv::cumulative_time, 0, 0};
}
apt::array<real_t, 3> Pdot_photon_pair_creation(pgn::Index_t I,
                                                const pgn::Grid_t& grid,
                                                const pgn::Field<3>&,
                                                const pgn::Field<3>&,
                                                const pgn::JField&) {
  const auto& pc = gv::pc_counter;
  const auto& dt = gv::cumulative_time;
  return {pc[1](I) / dt, pc[2](I) / dt, pc[3](I) / dt};
}

apt::array<real_t, 3> photon_emission_rate(pgn::Index_t I,
                                           const pgn::Grid_t& grid,
                                           const pgn::Field<3>&,
                                           const pgn::Field<3>&,
                                           const pgn::JField&) {
  auto x = gv::em_counter[0](I);
  return {x / gv::cumulative_time, 0, 0};
}
apt::array<real_t, 3> Pdot_photon_emission(pgn::Index_t I,
                                           const pgn::Grid_t& grid,
                                           const pgn::Field<3>&,
                                           const pgn::Field<3>&,
                                           const pgn::JField&) {
  const auto& em = gv::em_counter;
  const auto& dt = gv::cumulative_time;
  return {em[1](I) / dt, em[2](I) / dt, em[3](I) / dt};
}

template <particle::species SP>
apt::array<real_t, 3> J_by_species(pgn::Index_t I, const pgn::Grid_t& grid,
                                   const pgn::Field<3>&, const pgn::Field<3>&,
                                   const pgn::JField&) {
  // NOTE due to interpolation, the exported Jsp's don't sum up to J.
  const auto& Jsp = (*gv::Jsp)[SP];
  return msh::interpolate(Jsp, I2std(I), ShapeF());
}

apt::array<real_t, 3> skin_depth(pgn::Index_t I, const pgn::Grid_t& grid,
                                 const pgn::Field<3>&, const pgn::Field<3>&,
                                 const pgn::JField&) {
  return {msh::interpolate(*(gv::skin_depth), I2std(I), ShapeF())[0], 0, 0};
}

apt::array<real_t, 3> ptc_num(
    const particle::Properties& prop,
    const pgn::ParticleArray_t::const_particle_type& ptc) {
  return {1.0, 0.0, 0.0};
}

// NOTE energy should factor in mass
apt::array<real_t, 3> ptc_energy(
    const particle::Properties& prop,
    const pgn::ParticleArray_t::const_particle_type& ptc) {
  if (prop.mass_x > 0.01_r)
    return {prop.mass_x * std::sqrt(1.0_r + apt::sqabs(ptc.p())), 0.0, 0.0};
  else
    return {apt::abs(ptc.p()), 0.0, 0.0};
}

apt::array<real_t, 3> ptc_momentum(
    const particle::Properties& prop,
    const pgn::ParticleArray_t::const_particle_type& ptc) {
  if (prop.mass_x > 0.01_r)
    return {prop.mass_x * ptc.p(0), prop.mass_x * ptc.p(1),
            prop.mass_x * ptc.p(2)};
  else
    return {ptc.p(0), ptc.p(1), ptc.p(2)};
}

}  // namespace io

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
  gv::wpic2 = 2 * gv::Omega * gv::mu / conf["Np"].as<real_t>();

  const real_t gamma_fd = conf["pairs"]["gamma_fd"].as<real_t>();
  const real_t E_ph = conf["pairs"]["E_ph"].as<real_t>();
  const real_t gravity_strength = conf["forces"]["gravity"].as<real_t>();
  const real_t landau0_B_thr =
      conf["forces"]["landau0_ratio"].as<real_t>() * gv::mu;
  const real_t magnetic_convert_B_thr =
      conf["pairs"]["photon"]["magnetic_convert_ratio"].as<real_t>() * gv::mu;
  const std::array<real_t, 2> mfp = {
      conf["pairs"]["photon"]["mfp"][0].as<real_t>(),
      conf["pairs"]["photon"]["mfp"][1]
          .as<real_t>()};  // TODO this should use toml array parsing

  const real_t v_th = conf["atmosphere"]["v_th"].as<real_t>();
  const real_t atm_x = conf["atmosphere"]["multiplicity"].as<real_t>();

  const real_t damping_thickness = conf["damping"]["thickness"].as<real_t>();
  const real_t damping_rate = conf["damping"]["rate"].as<real_t>();
  gv::spinup_time = conf["spinup_time"].as<real_t>();

  const real_t r_e = gv::wpic2 * apt::dV(gv::supergrid) / (4.0 * M_PI);

  const int downsample_ratio =
      conf["schedules"]["export"]["downsample_ratio"].as_or<int>(1);

  const pgn::ExportSchedule sch_export = [&] {
    pgn::ExportSchedule res;
    const auto& cf = conf["schedules"]["export"];
    res.on = cf["on"].as_or<bool>(true);
    res.start = cf["start"].as_or<int>(0);
    res.interval, cf["interval"].as<int>();
    res.num_files = cf["num_files"].as_or<int>(1);
    return res;
  }();

  pgn::SimulationBuilder_t builder(args);

  builder.initialize_this_run_dir(datadir_prefix, project_name)
      .set_supergrid(gv::supergrid)
      .set_field_guard(gv::guard)
      .set_cartesian_topology(dims, periodic);

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

    builder.add_field_action<field::Axisymmetric>()
        .is_upper_axis(false)
        .set_name("AxisymmetrizeEBLower")
        .set_range(0, {0, gv::supergrid[0].dim()})
        .set_range(
            1, {-gv::guard, 1});  // NOTE 1 so as to set values right on axis

    builder.add_field_action<field::Axisymmetric>()
        .is_upper_axis(true)
        .set_name("AxisymmetrizeEBHigher")
        .set_range(0, {0, gv::supergrid[0].dim()})
        .set_range(
            1, {gv::supergrid[1].dim(), gv::supergrid[1].dim() + gv::guard});
  }

  // TODO
  {
    constexpr real_t gamma_off = 15.0;
    constexpr real_t Ndot_fd = 0.25;

    auto properties = set_up_particle_properties(
        gravity_strength, landau0_B_thr, gamma_fd, gamma_off, Ndot_fd, E_ph,
        magnetic_convert_B_thr, mfp);
    builder.set_particle_properties(properties);
  }

  {  // set up particle actions

    {
      int i_equator = gv::supergrid[1].dim() / 2 +
                      1;  // pick the cell whose lb is at equator.
      int half_width = (30.0_deg) / gv::supergrid[1].delta();
      builder.add_particle_action<particle::NoPhotonZone>()
          .set_interval(100)
          .set_name("NoPhotonZone")
          .set_range(0, {gv::supergrid[0].csba(std::log(18.0)),
                         gv::supergrid[0].dim() + gv::guard})
          .set_range(1, {i_equator - half_width, i_equator + half_width});

      builder.add_particle_action<particle::Annihilator>()
          .set_name("Annihilation")
          .set_range(0, {gv::supergrid[0].csba(std::log(18.0)),
                         gv::supergrid[0].dim() + gv::guard})
          .set_range(1, {i_equator - half_width, i_equator + half_width})
          .apply([&](auto& act) {
            auto& sch = act.schedule();
            const auto& cf = conf["schedules"]["annihilation"];
            sch.on = cf["on"].as_or<bool>(false);
            sch.start = cf["start"].as_or<int>(0);
            sch.interval = cf["interval"].as<int>();
          });
    }

    {
      const int escape_begin_index = gv::supergrid[0].csba(
          std::log(std::exp(gv::supergrid[0].upper()) -
                   0.5 * damping_thickness));  // at halfway so as to reduce
                                               // reflection of artifacts
      builder.add_particle_action<particle::Escaping>()
          .set_name("Escaping")
          .set_range(0,
                     {escape_begin_index, gv::supergrid[0].dim() + gv::guard})
          .set_range(1, {0, gv::supergrid[1].dim()});
    }

    builder.add_particle_action<particle::MainUpdater>()
        .set_name("MainUpdate")  // TODO set range??
        .set_update_q(
            Metric::geodesic_move<apt::vVec<real_t, 3>, apt::vVec<real_t, 3>>)
        .set_export_schedule(sch_export);

    {
      real_t N_atm_floor = std::exp(1.0) * 2.0 * gv::Omega * gv::mu / gv::wpic2;
      builder.add_particle_action<particle::Atmosphere>()
          .set_name("Atmosphere")
          .set_range(0, {gv::star_interior - 1, gv::star_interior})
          .set_range(1, {0, gv::supergrid[1].dim()})
          .set_thermal_velocity(v_th)
          .set_number_in_atmosphere(atm_x * N_atm_floor)
          .set_positive_charge(species::ion)
          .set_negative_charge(species::electron)
          .set_omega_t(omega_spinup)
          .set_normal_direction(0);
    }

    builder.add_particle_action<particle::Axisymmetric>()
        .set_name("AxisymmetrizeJLower")
        .set_range(0, {0, gv::supergrid[0].dim()})
        .set_range(1, {-gv::guard, 0})  // NOTE must not use 1 in place of 0
        .is_upper_axis(false);

    builder.add_particle_action<particle::Axisymmetric>()
        .set_name("AxisymmetrizeJHigher")
        .set_range(0, {0, gv::supergrid[0].dim()})
        .set_range(1,
                   {gv::supergrid[1].dim(), gv::supergrid[1].dim() + gv::guard})
        .is_upper_axis(true);

    builder.add_particle_action<particle::NewPtcAnalyzer>().set_name(
        "NewParticleAnalysis");  // TODO set range?
  }

  builder.add_init_cond_action<InitialCondition>()
      .set_range(0, {0, gv::supergrid[0].dim()})
      .set_range(1, {0, gv::supergrid[1].dim() +
                            1});  // NOTE +1 to include upper boundary

  builder.set_prior_export(io::do_prior_export)
      .set_post_export(io::do_post_export);

  auto add_common_exportees = [](auto& exporter) {
    exporter
        .add_exportee("E", 3, io::field_self<0>, io::average_when_downsampled)
        .add_exportee("B", 3, io::field_self<1>, io::average_when_downsampled)
        .add_exportee("J", 3, io::field_self<2>,
                      io::average_and_divide_flux_by_area)
        .add_exportee("PairCreationRate", 1, io::pair_creation_rate)
        .add_exportee("PdotPairCreation", 3, io::Pdot_photon_pair_creation)
        .add_exportee("PhotonEmissionRate", 1, io::photon_emission_rate)
        .add_exportee("PdotEmission", 3, io::Pdot_photon_emission);

    if (gv::skin_depth) {
      exporter.add_exportee("SkinDepth", 1, io::skin_depth,
                            io::average_when_downsampled);
    }

    if (gv::Jsp) {
      using namespace particle;
      for (auto sp : *gv::Jsp) {
        switch (sp) {
          case species::electron:
            exporter.add_exportee("Je", 3, io::J_by_species<species::electron>,
                                  io::average_and_divide_flux_by_area);
            break;
          case species::positron:
            exporter.add_exportee("Jp", 3, io::J_by_species<species::positron>,
                                  io::average_and_divide_flux_by_area);
            break;
          case species::ion:
            exporter.add_exportee("Ji", 3, io::J_by_species<species::ion>,
                                  io::average_and_divide_flux_by_area);
            break;
          default:;
        }
      }
    }
    exporter.add_exportee("Num", 1, io::ptc_num)
        .add_exportee("E", 1, io::ptc_energy)
        .add_exportee("P", 3, io::ptc_momentum);
  };

  builder.add_exporter()
      .set_is_collinear_mesh(false)
      .set_downsample_ratio(downsample_ratio)
      .set_data_dir("")
      .set_range({{{0, gv::supergrid[0].dim() + gv::guard},
                   {0, gv::supergrid[1].dim() + 1}}})
      .apply(add_common_exportees);

  {
    int i_equator = gv::supergrid[1].dim() / 2 +
                    1;  // pick the cell whose lb is at equator.
    int half_width = (30.0_deg) / gv::supergrid[1].delta();
    builder.add_exporter()
        .set_is_collinear_mesh(false)
        .set_downsample_ratio(1)
        .set_data_dir("current_sheet")
        .set_range({{{gv::supergrid[0].csba(std::log(4.5)),
                      gv::supergrid[0].dim() + gv::guard},
                     {i_equator - half_width, i_equator + half_width}}})
        .apply(add_common_exportees);
  }

  builder.add_exporter()
      .set_is_collinear_mesh(false)
      .set_downsample_ratio(1)
      .set_data_dir("Br_full")
      .set_range({{{0, gv::supergrid[0].dim() + gv::guard},
                   {0, gv::supergrid[1].dim() + 1}}})
      .add_exportee("Br", 1, io::Br_alone);

  builder.set_print_timestep_to_stdout_interval(
      conf["schedules"]["print_timestep_to_stdout_interval"].as_or<int>(100));

  builder.set_scattering_data_in_vitals(&gv::N_scat);

  {
    auto& sch = builder.sort_ptcs_schedule();
    const auto& cf = conf["schedules"]["sort"];
    sch.on = cf["on"].as_or<bool>(true);
    sch.start = cf["start"].as_or<int>(0);
    sch.interval = cf["interval"].as<int>();
  }

  {
    auto& sch = builder.export_schedule();
    sch = sch_export;
  }

  {
    auto& sch = builder.checkpoint_schedule();
    const auto& cf = conf["schedules"]["checkpoint"];
    sch.on = cf["on"].as_or<bool>(false);
    sch.start = cf["start"].as_or<int>(1);
    sch.interval = cf["interval"].as<int>();
    sch.num_files = cf["num_files"].as_or<int>(1);
    sch.max_num_checkpoints = cf["max_num_checkpoints"].as_or<int>(1);
  }

  {
    auto& sch = builder.load_balancing_schedule();
    const auto& cf = conf["schedules"]["load_balance"];
    sch.on = cf["on"].as_or<bool>(true);
    sch.start = cf["start"].as_or<int>(0);
    sch.interval = cf["interval"].as<int>();
    sch.target_load = cf["target_load"].as_or<int>(100000);
  }

  {
    auto& sch = builder.profiling_schedule();
    const auto& cf = conf["schedules"]["profiling"];
    sch.on = cf["on"].as_or<bool>(false);
    sch.start = cf["start"].as_or<int>(0);
    sch.interval = cf["interval"].as<int>();
    sch.max_entries = cf["max_entries"].optional<int>();
  }

  {
    auto& sch = builder.vitals_schedule();
    const auto& cf = conf["schedules"]["vitals"];
    sch.on = cf["on"].as_or<bool>(true);
    sch.start = cf["start"].as_or<int>(0);
    sch.interval = cf["interval"].as<int>();
  }

  builder.set_extra_init([](const auto& properties, const auto& localgrid) {
    for (auto sp : properties) gv::N_scat.insert(sp, 0);

    if (false) {  // activate Jsp
      gv::Jsp.emplace();
      for (auto sp : properties) (*gv::Jsp).insert(sp, {});
    }

    if (false) {  // activate skin_depth
      gv::skin_depth.emplace();
    }

    pgn::Index_t bulk_dims;
    for (int i = 0; i < DGrid; ++i) bulk_dims[i] = localgrid[i].dim();
    // FIXME range with no guard gives memory error. Interpolation in export
    // needs them.
    auto range = apt::make_range({}, bulk_dims, gv::guard);
    gv::pc_counter = {range};
    gv::em_counter = {range};
  });

  auto& sim = builder.build();

  // TODO having to have user call this is a bit error_prone
  const auto init_ts = sim.initial_timestep();
  if (mpi::world.rank() == 0) std::cout << "Launch" << std::endl;
  for (int ts = init_ts; ts < init_ts + n_timesteps; ++ts) {
    sim.evolve(ts, dt);
  }

  return 0;
}
