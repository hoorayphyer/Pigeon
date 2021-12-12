#pragma once
#include <functional>
#include <memory>
#include <type_traits>
#include <vector>

#include "io/data_exporter.hpp"
#include "simulator/action.hpp"
#include "simulator/simulator.hpp"

namespace pic {

struct CLIArgs;

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
class SimulationBuilder {
 public:
  using FieldAction_t = FieldAction<DGrid, R, RJ>;
  using ParticleAction_t = ParticleAction<DGrid, R, S, RJ>;
  using InitialConditionAction_t = InitialConditionAction<DGrid, R, S, RJ>;
  using PostResumeAction_t = PostResumeAction<DGrid, R, S, RJ>;
  using DataExporter_t = io::DataExporter<RD, DGrid, R, S, RJ>;
  using ExportBundle_t = ExportBundle<DGrid, R, S, RJ>;
  using Simulator_t = Simulator<DGrid, R, S, RJ, RD>;

  SimulationBuilder(CLIArgs args);
  ~SimulationBuilder();

  /**
     add a field action of ConcreteAction constructed with `ctor_args`.

     @return reference to the newly added action of type ConcreteAction
   */
  template <typename ConcreteAction, typename... Args>
  ConcreteAction& add_field_action(Args&&... ctor_args) {
    return add_action_impl<FieldAction_t, ConcreteAction>(
        std::forward<Args>(ctor_args)...);
  }

  /**
     add a particle action of ConcreteAction constructed with `ctor_args`.

     @return reference to the newly added action of type ConcreteAction
   */
  template <typename ConcreteAction, typename... Args>
  ConcreteAction& add_particle_action(Args&&... ctor_args) {
    return add_action_impl<ParticleAction_t, ConcreteAction>(
        std::forward<Args>(ctor_args)...);
  }

  /**
     add an action of ConcreteAction constructed with `ctor_args` for setting up
     initial condition actions

     @return reference to the newly added action of type ConcreteAction
   */
  template <typename ConcreteAction, typename... Args>
  ConcreteAction& add_init_cond_action(Args&&... ctor_args) {
    return add_action_impl<InitialConditionAction_t, ConcreteAction>(
        std::forward<Args>(ctor_args)...);
  }

  /**
     add an action of ConcreteAction constructed with `ctor_args` for setting up
     post resume actions

     @return reference to the newly added action of type ConcreteAction
   */
  template <typename ConcreteAction, typename... Args>
  ConcreteAction& add_post_resume_action(Args&&... ctor_args) {
    return add_action_impl<PostResumeAction_t, ConcreteAction>(
        std::forward<Args>(ctor_args)...);
  }

  /**
     add an exporter.

     @return reference to the newly added exporter
   */
  DataExporter_t& add_exporter() { return m_sim->m_exporters.emplace_back(); }

  SimulationBuilder& set_extra_init(
      std::function<void(const particle::map<particle::Properties>& properties,
                         const apt::Grid<R, DGrid>& localgrid)>
          f) {
    m_sim->m_f_extra_init = std::move(f);
    return *this;
  }

  SimulationBuilder& set_prior_export(
      std::function<void(const ExportBundle_t&)> f) {
    m_sim->m_f_prior_export.emplace(std::move(f));
    return *this;
  }

  SimulationBuilder& set_post_export(
      std::function<void(const ExportBundle_t&)> f) {
    m_sim->m_f_post_export.emplace(std::move(f));
    return *this;
  }

  SimulationBuilder& set_scattering_data_in_vitals( const particle::map<R>* N_scat) {
    m_sim->m_N_scat = N_scat;
    return *this;
  }

  SimulationBuilder& initialize_this_run_dir(std::string prefix,
                                             std::string project_name);

  SimulationBuilder& set_supergrid(const apt::Grid<R, DGrid>& supergrid) {
    m_supergrid.emplace(supergrid);
    return *this;
  }

  SimulationBuilder& create_cartesian_topology(
      const apt::array<int, DGrid>& dims,
      const apt::array<bool, DGrid>& periodic);

  SimulationBuilder& set_particle_properties(
      particle::map<particle::Properties> props) {
    m_props = std::move(props);
    return *this;
  }

  SimulationBuilder& set_field_guard(int guard) {
    m_fld_guard = guard;
    return *this;
  }

  auto& sort_ptcs_schedule() { return m_sim->m_sch_sort_ptcs; }

  auto& export_schedule() { return m_sim->m_sch_export; }

  auto& checkpoint_schedule() { return m_sim->m_sch_ckpt; }

  auto& load_balancing_schedule() { return m_sim->m_sch_dlb; }

  auto& profiling_schedule() { return m_sim->m_sch_prof; }

  auto& vitals_schedule() { return m_sim->m_sch_vitals; }

  SimulationBuilder& set_print_timestep_to_stdout_interval(int interval) {
    m_sim->m_print_timestep_to_stdout_interval = interval;
    return *this;
  }

  Simulator_t& build();

  // TODO 1. maybe use concept 2. can we avoid having to have user call this
  template <typename Ptc>
  SimulationBuilder& commit_particle_type_for_mpi() {
    mpi::commit(mpi::Datatype<Ptc>{});
    m_is_mpi_particle_committed = true;
    return *this;
  }

 private:
  CLIArgs m_args;
  bool m_is_mpi_particle_committed = false;
  bool m_is_sim_complete = false;

  // these are with global ranges;
  std::vector<std::unique_ptr<InitialConditionAction_t>> m_ic_actions;
  std::vector<std::unique_ptr<PostResumeAction_t>> m_post_resume_actions;

  std::optional<std::string> m_this_run_dir;

  std::optional<apt::Grid<R, DGrid>> m_supergrid;
  std::optional<std::optional<mpi::CartComm>> m_cart;
  apt::array<int, DGrid> m_dims;
  apt::array<bool, DGrid> m_periodic;
  particle::map<particle::Properties> m_props;

  std::optional<int> m_total_timesteps;
  std::optional<int> m_fld_guard;

  // use std::optional so as to destruct m_sim before mpi::finalize, see
  // destructor.
  std::optional<Simulator_t> m_sim;

  template <typename Action>
  auto& get_actions() {
    if constexpr (std::is_same_v<Action, FieldAction_t>) {
      return m_sim->m_fld_actions;
    } else if constexpr (std::is_same_v<Action, ParticleAction_t>) {
      return m_sim->m_ptc_actions;
    } else if constexpr (std::is_same_v<Action, InitialConditionAction_t>) {
      return m_ic_actions;
    } else if constexpr (std::is_same_v<Action, PostResumeAction_t>) {
      return m_post_resume_actions;
    } else {
      []<bool flag = false>() { static_assert(flag, "Unknown action type!"); }
      ();
    }
  }

  template <typename Action, typename ConcreteAction, typename... Args>
  ConcreteAction& add_action_impl(Args&&... ctor_args) {
    static_assert(std::is_base_of_v<Action, ConcreteAction>);
    auto& actions = get_actions<Action>();
    auto& uniq_ptr = actions.emplace_back(
        new ConcreteAction(std::forward<Args>(ctor_args)...));
    return static_cast<ConcreteAction&>(*uniq_ptr);
  }

  std::string precondition() const;

  int load_init_cond(Simulator_t& sim);

  void set_up_journal() const;
  void populate_this_run_dir() const;
};
}  // namespace pic
