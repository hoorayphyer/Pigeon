#include <algorithm>
#include <stdexcept>

#include "simulator/builder.hpp"

namespace pic {

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
SimulationBuilder<DGrid, R, S, RJ, RD>&
SimulationBuilder<DGrid, R, S, RJ, RD>::create_cartesian_topology(
    const apt::array<int, DGrid>& dims,
    const apt::array<bool, DGrid>& periodic) {
  m_cart.emplace();
  const int nprmy =
      std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<int>());

  if (nprmy > mpi::world.size())
    throw std::invalid_argument(
        "Size of the Cartesian topology exceeds the size of world!");

  // simply use first few members in mpi::world as primaries
  bool is_prmy = mpi::world.rank() < nprmy;
  auto comm_tmp = mpi::world.split({is_prmy}, mpi::world.rank());
  if (is_prmy) {
    std::vector<int> d(DGrid);
    std::vector<bool> p(DGrid);
    for (int i = 0; i < DGrid; ++i) {
      d[i] = dims[i];
      p[i] = periodic[i];
    }
    m_cart->emplace(*comm_tmp, d, p);
  }

  m_dims = dims;
  m_periodic = periodic;

  return *this;
}

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
std::string SimulationBuilder<DGrid, R, S, RJ, RD>::precondition() const {
  std::string msg = "";
  if (!m_supergrid) {
    msg += "- must set supergrid\n";
  }
  if (!m_cart) {
    msg += "- must initialize cartesian topology\n";
  }
  if (!m_this_run_dir) {
    msg += "- must set this_run_dir\n";
  }
  if (m_ic_actions.size() == 0) {
    msg += "- no initial condition is given\n";
  }
  if (!m_total_timesteps) {
    msg += "- must set the total number of timesteps\n";
  }
  if (!m_fld_guard) {
    msg += "- must set the field guard\n";
  }
  return msg;
}

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
int SimulationBuilder<DGrid, R, S, RJ, RD>::load_init_cond(Simulator_t& sim) {
  int init_ts = 0;
  // TODO profiling_plan
  // if (profiling_plan.on and profiling_plan.is_qualified()) {
  //   lgr::file << "--- Loading Initial Condition..." << std::endl;
  // }

  if (m_checkpoint_dir) {
    init_ts = ckpt::load_checkpoint(*m_checkpoint_dir, sim.m_ens_opt,
                                    sim.m_cart_opt, sim.m_E, sim.m_B,
                                    sim.m_particles, sim.m_properties);
    if (sim.m_ens_opt) {
      sim.update_parts();

      auto bundle = typename PostResumeAction_t::Bundle_t{
          sim.m_E,    sim.m_B,       sim.m_J, sim.m_particles, sim.m_properties,
          sim.m_grid, sim.m_ens_opt, init_ts, *m_this_run_dir};

      for (auto& act : m_post_resume_actions) {
        // TODO also check range emptiness
        if (!act) continue;
        sim.taylor(act->ranges());
        (*act)(bundle);
      }
    }
    ++init_ts;  // checkpoint is saved at the end of a timestep
  } else {
    // FIXME a temporary fix, which may crash under the edge case in which
    // initially many particles are created.
    if (sim.m_cart_opt) {
      auto bundle = typename InitialConditionAction_t::Bundle_t{
          sim.m_E,         sim.m_B,          sim.m_J,
          sim.m_particles, sim.m_properties, sim.m_grid};
      for (auto& act : m_ic_actions) {
        // TODO also check range emptiness
        if (!act) continue;
        sim.taylor(act->ranges());
        (*act)(bundle);
      }

      for (auto sp : sim.m_particles) {
        for (auto ptc : sim.m_particles[sp]) {  // TODOL semantics
          ptc.set(sp);
        }
      }
      field::copy_sync_guard_cells(sim.m_E, *sim.m_cart_opt);
      field::copy_sync_guard_cells(sim.m_B, *sim.m_cart_opt);
      field::copy_sync_guard_cells(sim.m_J, *sim.m_cart_opt);
    }
  }
  // broadcast to ensure uniformity
  mpi::world.broadcast(0, &init_ts, 1);

  // TODO the following can be made into a separate action
  // // initialize trace_counters to ensure unique trace serial numbers across
  // // runs
  // pic::Traman::init(mpi::world.rank(), _particles);

  // TODO profiling plan
  // if (profiling_plan.on and profiling_plan.is_qualified()) {
  //   lgr::file << "--- Done." << std::endl;
  // }

  return init_ts;
}

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
SimulationBuilder<DGrid, R, S, RJ, RD>::Simulator_t
SimulationBuilder<DGrid, R, S, RJ, RD>::build() {
  if (m_is_build_called) {
    throw std::runtime_error(
        "must call build() no more than once on a SimulationBuilder object!");
  }

  {
    auto msg = precondition();
    if (!msg.empty()) {
      msg = "Error: incomplete simulation setup\n" + msg;
      throw std::runtime_error(msg);
    }
  }

  Simulator_t sim;
  sim.initialize(std::move(*m_supergrid), std::move(*m_cart),
                 std::move(m_props), std::move(m_periodic), m_dims);
  auto init_ts = load_init_cond(sim);

  sim.m_rng.set_seed(init_ts + mpi::world.rank());

  // TODO field/ptc actions and their orig ranges
  sim.m_f_prior_export = std::move(m_f_prior_export);
  sim.m_exporters = std::move(m_exporters);
  sim.m_f_post_export = std::move(m_f_post_export);
  sim.m_f_custom_step = std::move(m_f_custom_step);
  sim.m_fld_guard = *m_fld_guard;
  sim.m_this_run_dir = *m_this_run_dir;

  m_is_build_called = true;

  return sim;
}
}  // namespace pic