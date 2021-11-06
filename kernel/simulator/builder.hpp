#pragma once
#include <functional>
#include <memory>
#include <type_traits>
#include <vector>

#include "io/data_exporter.hpp"
#include "simulator/action.hpp"
#include "simulator/simulator.hpp"

namespace pic {

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
struct SimulationBuilder {
 public:
  using FieldAction_t = FieldAction<DGrid, R, RJ>;
  using ParticleAction_t = ParticleAction<DGrid, R, S, RJ>;
  using InitialConditionAction_t = InitialConditionAction<DGrid, R, S, RJ>;
  using PostResumeAction_t = PostResumeAction<DGrid, R, S, RJ>;
  using DataExporter_t = io::DataExporter<RD, DGrid, R, S, RJ>;
  using ExportBundle_t = ExportBundle<DGrid, R, S, RJ>;
  using Simulator_t = Simulator<DGrid, R, S, RJ, RD>;

  /**
     add a field action of ConcreteAction constructed with `ctor_args`.

     @return reference to the newly added action of type ConcreteAction
   */
  template <typename ConcreteAction, typename... Args>
  ConcreteAction& add_field_action(Args&&... ctor_args) {
    return add_action<FieldAction_t, ConcreteAction>(
        std::forward<Args>(ctor_args)...);
  }

  /**
     add a particle action of ConcreteAction constructed with `ctor_args`.

     @return reference to the newly added action of type ConcreteAction
   */
  template <typename ConcreteAction, typename... Args>
  ConcreteAction& add_particle_action(Args&&... ctor_args) {
    return add_action<ParticleAction_t, ConcreteAction>(
        std::forward<Args>(ctor_args)...);
  }

  /**
     add an action of ConcreteAction constructed with `ctor_args` for setting up
     initial condition actions

     @return reference to the newly added action of type ConcreteAction
   */
  template <typename ConcreteAction, typename... Args>
  ConcreteAction& add_init_cond_action(Args&&... ctor_args) {
    return add_action<InitialConditionAction_t, ConcreteAction>(
        std::forward<Args>(ctor_args)...);
  }

  /**
     add an action of ConcreteAction constructed with `ctor_args` for setting up
     post resume actions

     @return reference to the newly added action of type ConcreteAction
   */
  template <typename ConcreteAction, typename... Args>
  ConcreteAction& add_post_resume_action(Args&&... ctor_args) {
    return add_action<PostResumeAction_t, ConcreteAction>(
        std::forward<Args>(ctor_args)...);
  }

  /**
     add an exporter.

     @return reference to the newly added exporter
   */
  DataExporter_t& add_exporter() { return m_exporters.emplace_back(); }

  auto& add_extra_init(std::function<void()> f) {
    m_f_extra_init = std::move(f);
    return *this;
  }

  auto& set_prior_export(std::function<void(const ExportBundle_t&)> f) {
    m_f_prior_export = std::move(f);
    return *this;
  }

  auto& set_post_export(std::function<void(const ExportBundle_t&)> f) {
    m_f_post_export = std::move(f);
    return *this;
  }

  auto& add_custom_step(std::function<void()> f) {
    m_f_custom_step = std::move(f);
    return *this;
  }

  auto& set_this_run_dir(std::string dir) {
    m_this_run_dir = std::move(dir);
    return *this;
  }

  // TODO
  // add_supergrid();
  // add_ptc_prop();
  // init_cart. See pic.cpp, make_cart
  // set_checkpoint_dir();
  // set_periodic();

  Simulator_t build();

 private:
  // these are with global ranges;
  std::vector<std::unique_ptr<FieldAction_t>> m_fld_actions;
  std::vector<std::unique_ptr<ParticleAction_t>> m_ptc_actions;
  std::vector<std::unique_ptr<InitialConditionAction_t>> m_ic_actions;
  std::vector<std::unique_ptr<PostResumeAction_t>> m_post_resume_actions;
  std::vector<DataExporter_t> m_exporters;

  std::optional<std::function<void()>> m_f_extra_init;
  std::optional<std::function<void(const ExportBundle_t&)>> m_f_prior_export;
  std::optional<std::function<void(const ExportBundle_t&)>> m_f_post_export;
  std::optional<std::function<void()>> m_f_custom_step;

  std::optional<std::string> m_this_run_dir;

  std::optional<apt::Grid<R, DGrid>> m_supergrid;
  std::optional<std::optional<mpi::CartComm>> m_cart;
  particle::map<particle::Properties> m_props;

  std::optional<std::string> m_checkpoint_dir;
  std::optional<int> m_total_timesteps;
  std::optional<apt::array<bool, DGrid>> m_periodic;

  template <typename Action>
  auto& get_actions() {
    if constexpr (std::is_same_v<Action, FieldAction_t>) {
      return m_fld_actions;
    } else if constexpr (std::is_same_v<Action, ParticleAction_t>) {
      return m_ptc_actions;
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
};
}  // namespace pic
