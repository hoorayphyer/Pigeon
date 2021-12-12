#pragma once

#include <memory>
#include <unordered_map>

#include "ckpt/autosave.hpp"
#include "ckpt/checkpoint.hpp"
#include "dye/dynamic_balance.hpp"
#include "field/sync.hpp"
#include "field/yee.hpp"
#include "filesys/filesys.hpp"
#include "logger/logger.hpp"
#include "particle/sorter.hpp"
#include "timer/timer.hpp"

#if PIC_DEBUG
#include "debug/debugger.hpp"
#endif

#include "io/data_exporter.hpp"
#include "simulator/action.hpp"
#include "simulator/schedule.hpp"

namespace pic {

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
class SimulationBuilder;

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
struct Simulator {
 private:
  using FieldAction_t = FieldAction<DGrid, R, RJ>;
  using ParticleAction_t = ParticleAction<DGrid, R, S, RJ>;
  using DataExporter_t = io::DataExporter<RD, DGrid, R, S, RJ>;
  using ExportBundle_t = ExportBundle<DGrid, R, S, RJ>;

  apt::Grid<R, DGrid> m_supergrid;
  std::optional<mpi::CartComm> m_cart_opt;
  util::Rng<R> m_rng;
  int m_fld_guard;
  std::string m_this_run_dir;
  int m_initial_timestep = 0;

  apt::Grid<R, DGrid> m_grid;
  std::optional<dye::Ensemble<DGrid>> m_ens_opt;

  particle::map<particle::Properties> m_properties;
  field::Field<R, 3, DGrid> m_E;
  field::Field<R, 3, DGrid> m_B;
  field::Field<RJ, 3, DGrid> m_J;
  particle::map<particle::array<R, S>> m_particles;

  std::vector<apt::array<apt::Range, DGrid>> m_fld_action_orig_ranges;
  std::vector<std::unique_ptr<FieldAction_t>> m_fld_actions;
  std::vector<apt::array<apt::Range, DGrid>> m_ptc_action_orig_ranges;
  std::vector<std::unique_ptr<ParticleAction_t>> m_ptc_actions;

  std::vector<particle::Particle<R, S>> m_ptc_buffer;
  apt::array<bool, DGrid> m_periodic;

  std::optional<std::function<void(const ExportBundle_t&)>> m_f_prior_export;
  std::vector<apt::array<apt::Range, DGrid>> m_exporters_orig_ranges;
  std::vector<DataExporter_t> m_exporters;
  std::optional<std::function<void(const ExportBundle_t&)>> m_f_post_export;

  std::optional<
      std::function<void(const particle::map<particle::Properties>& properties,
                         const apt::Grid<R, DGrid>& localgrid)>>
      m_f_extra_init;
  const particle::map<R>* m_N_scat = nullptr;

  Schedule m_sch_sort_ptcs;
  ExportSchedule m_sch_export;
  CheckpointSchedule m_sch_ckpt;
  DynamicLoadBalanceSchedule m_sch_dlb;
  ProfilingSchedule m_sch_prof;
  int m_print_timestep_to_stdout_interval = 100;
  Schedule m_sch_vitals;

  Simulator() = default;

  void taylor(apt::array<apt::Range, DGrid>& a) const;

  void update_parts();

  void print_vitals(const std::string& filename, R t_phys) const;

 public:
  friend class SimulationBuilder<DGrid, R, S, RJ, RD>;

  int initial_timestep() const noexcept { return m_initial_timestep; }

  void evolve(int timestep, R dt);
};
}  // namespace pic
