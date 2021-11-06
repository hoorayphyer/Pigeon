#include <stdexcept>

#include "simulator/builder.hpp"

namespace pic {
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
    msg += "- no initial condition is given\n"
  }
  if (!m_total_timesteps) {
    msg += "- must set the total number of timesteps\n";
  }
  if (!m_periodic) {
    msg += "- must set periodicity\n";
  }
  return msg;
}

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
int SimulationBuilder<DGrid, R, S, RJ, RD>::load_init_cond(Simulator_t& sim) {
  int init_ts = 0;
  if (profiling_plan.on and profiling_plan.is_qualified()) {
    lgr::file << "--- Loading Initial Condition..." << std::endl;
  }

  if (m_checkpoint_dir) {
    init_ts = ckpt::load_checkpoint(*m_checkpoint_dir, sim.m_ens_opt,
                                    sim.m_cart_opt, sim.m_E, sim.m_B,
                                    sim.m_particles, sim.m_properties);
    if (sim.m_ens_opt) {
      sim.update_parts(*m_ens_opt);

      PostResumeAction_t::Bundle_t bundle{
          sim.m_E,    sim.m_B,       sim.m_J, sim.m_particles, sim.m_properties,
          sim.m_grid, sim.m_ens_opt, init_ts, *m_this_run_dir};

      for (auto& act : m_post_resume_actions) {
        sim.taylor(act, *m_periodic);
        act(bundle);
      }
    }
    ++init_ts;  // checkpoint is saved at the end of a timestep
  } else {
    // FIXME a temporary fix, which may crash under the edge case in which
    // initially many particles are created.
    if (sim.m_cart_opt) {
      InitialConditionAction_t::Bundle_t bundle{
          sim.m_E,         sim.m_B,          sim.m_J,
          sim.m_particles, sim.m_properties, sim.m_grid};
      for (auto& act : m_ic_actions) {
        sim.taylor(act, *m_periodic);
        act(bundle);
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

  // TODO
  // // initialize trace_counters to ensure unique trace serial numbers across
  // // runs
  // pic::Traman::init(mpi::world.rank(), _particles);

  if (profiling_plan.on and profiling_plan.is_qualified()) {
    lgr::file << "--- Done." << std::endl;
  }

  return init_ts;
}

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
SimulationBuilder<DGrid, R, S, RJ, RD>::Simulator_t build() {
  {
    auto msg = precondition();
    if (!msg.empty()) {
      msg = "Error: incomplete simulation setup\n" + msg;
      throw std::runtime_error(msg);
    }
  }

  Simulator_t sim;
  sim.initialize(std::move(*m_supergrid), std::move(*m_cart),
                 std::move(m_props));
  auto init_ts = load_init_cond(sim);

  sim.m_rng.set_seed(init_ts + mpi::world.rank());

  return sim;
}
}  // namespace pic
