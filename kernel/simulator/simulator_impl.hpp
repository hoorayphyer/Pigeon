#include <cassert>

#include "simulator/simulator.hpp"

namespace pic {

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
void Simulator<DGrid, R, S, RJ, RD>::taylor(
    apt::array<apt::Range, DGrid>& a) const {
  // NOTE range is assumed to be noempty [,)
  auto f = [](int Ib_global, int Ie_global, const apt::Grid1D<R>& supergrid,
              const apt::Grid1D<R>& localgrid,
              bool is_periodic) noexcept -> apt::pair<int> {
    // shift is the index of local lower in the global grid
    int shift = static_cast<int>(
        (localgrid.lower() - supergrid.lower()) / localgrid.delta() + 0.5);
    Ib_global -= shift;
    Ie_global -= shift;
    // now Ib_ and Ie_global are with respect to the current local grid

    // Extend to infinity on true boundaries
    // NOTE there may be actions done completely in the guard cells, such as
    // assigining values
    int lb_local = 0;
    int ub_local = localgrid.dim();
    if (!is_periodic) {
      if (std::abs(localgrid.lower() - supergrid.lower()) < supergrid.delta())
        lb_local = std::numeric_limits<int>::min();
      if (std::abs(localgrid.upper() - supergrid.upper()) < supergrid.delta())
        ub_local = std::numeric_limits<int>::max();
    }

    if (Ie_global <= lb_local || Ib_global >= ub_local) {
      // range not applicable on current local patch
      return {0, 0};
    } else {
      Ib_global = std::max<int>(lb_local, Ib_global);
      Ie_global = std::min<int>(Ie_global, ub_local);
      return {Ib_global, Ie_global};
    }
  };

  for (int i = 0; i < DGrid; ++i) {
    auto [b_new, e_new] =
        f(a[i].begin(), a[i].end(), m_supergrid[i], m_grid[i], m_periodic[i]);
    a[i].begin() = b_new;
    a[i].end() = e_new;
  }
}

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
void Simulator<DGrid, R, S, RJ, RD>::update_parts() {
  assert(m_ens_opt);
  const auto& ens = *m_ens_opt;
  for (int i = 0; i < DGrid; ++i) {
    m_grid[i] =
        m_supergrid[i].divide(ens.cart_topos[i].dim(), ens.cart_coords[i]);
  }

  auto retaylor_ranges = [this](const auto& ranges_orig, auto& actions) {
    assert(actions.size() == ranges_orig.size());
    for (auto i = 0u; i < actions.size(); ++i) {
      if (!actions[i]) continue;
      auto& ranges = actions[i]->ranges();
      ranges = ranges_orig[i];
      taylor(ranges);
    }
  };

  retaylor_ranges(m_fld_action_orig_ranges, m_fld_actions);
  retaylor_ranges(m_ptc_action_orig_ranges, m_ptc_actions);
  {  // rataylor ranges for exporters
    const auto& ranges_orig = m_exporters_orig_ranges;
    auto& actions = m_exporters;
    assert(actions.size() == ranges_orig.size());
    for (auto i = 0u; i < actions.size(); ++i) {
      auto& ranges = actions[i].get_range();
      ranges = ranges_orig[i];
      taylor(ranges);
    }
  }

  if (m_f_extra_init) {
    (*m_f_extra_init)(m_properties, m_grid);
  }
}

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
void Simulator<DGrid, R, S, RJ, RD>::print_vitals(const std::string& filename,
                                                  R t_phys) const {
  static int counter = 0;
  constexpr int interval = 40;  // how often to print the table description row
  static tmr::Timestamp stopwatch;

  const auto& ens = *m_ens_opt;

  std::vector<double> buffer;
  {
    for (auto sp : m_particles) {
      double num = 0.0;
      for (const auto& ptc : m_particles[sp]) {
        if (ptc.is(particle::flag::exist)) num += ptc.frac();
      }
      buffer.push_back(num);
    }
    if (m_N_scat) {
      for (auto sp : *m_N_scat) buffer.push_back((*m_N_scat)[sp]);
    }
    ens.reduce_to_chief(mpi::by::SUM, buffer.data(), buffer.size());
  }
  if (!m_cart_opt) return;

  buffer.push_back(ens.size());
  m_cart_opt->template reduce<true>(mpi::by::SUM, 0, buffer.data(),
                                    buffer.size());

  if (m_cart_opt->rank() != 0) return;

  // significant only on world rank 0
  static double t_phys_prev = -1e-6;
  static std::vector<double> num_ptcs_prev(m_particles.size(), 0);
  static std::optional<std::vector<double>> num_scat_prev = [this] {
    std::optional<std::vector<double>> res;
    if (m_N_scat) {
      res.emplace(std::vector<double>((*m_N_scat).size(), 0));
    }
    return res;
  }();

  {
    std::ofstream out(filename, std::ios_base::app);
    if (counter % interval == 0) {
      out << "t_phys|\tnprocs|\tlapse/hr|\tTotal load(rate)";
      if (m_N_scat) {
        out << "|\tscattering creation rate";
      }
      out << std::endl;
      out << "species ordering : ";
      for (auto sp : m_particles) out << m_properties[sp].nickname << " ";
      if (m_N_scat) {
        // print separately in case it differs from above
        out << "|\t";
        for (auto sp : *m_N_scat) out << m_properties[sp].nickname << " ";
      }
      out << std::endl;
      counter = 0;
    }
    ++counter;
    out << apt::fmt("%.2f", t_phys) << "|\t" << buffer.back() << "|\t"
        << apt::fmt("%8.2f", stopwatch.lapse().in_units_of("s").val() / 3600.0)
        << "|\t";
    auto* p1 = buffer.data();
    for (int i = 0; i < m_particles.size(); ++i) {
      out << apt::fmt("%.2e", p1[i]) << "("
          << apt::fmt("%.2e",
                      (p1[i] - num_ptcs_prev[i]) / (t_phys - t_phys_prev))
          << ") ";
      num_ptcs_prev[i] = p1[i];
    }
    if (m_N_scat) {
      out << "|\t";
      auto* p2 = buffer.data() + m_particles.size();
      for (int i = 0; i < m_N_scat->size(); ++i) {
        out << apt::fmt("%.2e",
                        (p2[i] - (*num_scat_prev)[i]) / (t_phys - t_phys_prev))
            << " ";
        (*num_scat_prev)[i] = p2[i];
      }
    }
    out << std::endl;
    out.close();

    t_phys_prev = t_phys;
  }

  return;
}

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
void Simulator<DGrid, R, S, RJ, RD>::evolve(int timestep, R dt) {
  // NOTE we choose to do particle update before field update so that in
  // saving snapshots we don't need to save _J

  if (timestep % m_print_timestep_to_stdout_interval == 0 &&
      mpi::world.rank() == 0)
    std::cout << "==== Timestep " << timestep << " ====" << std::endl;

#if PIC_DEBUG
  debug::timestep = timestep;
  debug::world_rank = mpi::world.rank();
  if (_ens_opt) debug::ens_label = _ens_opt->label();
#endif

  std::optional<tmr::Timestamp> stamp_all;
  std::optional<tmr::Timestamp> stamp;
  if (m_sch_prof.is_do(timestep) and m_sch_prof.is_qualified()) {
    stamp_all.emplace();
    stamp.emplace();
    if (m_sch_prof.is_reached_max_entries(timestep)) lgr::file.clear();
    lgr::file << "==== Timestep " << timestep << " ====" << std::endl;
    lgr::file.indent_append("\t");
    if (m_ens_opt)
      lgr::file % "ensemble label = " << m_ens_opt->label() << std::endl;
    else
      lgr::file % "ensemble label = NONE" << std::endl;
  }

  auto TIMING = [&stamp](std::string message, auto f) {
    if (stamp) {
      lgr::file % message << "==>>" << std::endl;
      stamp.emplace();
    }
    f();
    if (stamp) {
      lgr::file % "\tLapse " << stamp->lapse().in_units_of("ms") << std::endl;
    }
  };
#define START [&]()

  if (m_ens_opt) {
    const auto& ens = *m_ens_opt;

    TIMING(
        "BroadcastEB", START {
          // TODOL reduce number of communications?
          for (int i = 0; i < 3; ++i)
            ens.intra.broadcast(ens.chief, m_E[i].data().data(),
                                m_E[i].data().size());
          for (int i = 0; i < 3; ++i)
            ens.intra.broadcast(ens.chief, m_B[i].data().data(),
                                m_B[i].data().size());
        });

    m_J.reset();

    if (m_sch_sort_ptcs.is_do(timestep)) {
      TIMING(
          "SortParticles", START {
            for (auto sp : m_particles) particle::sort(m_particles[sp], m_grid);
          });
    }

    if (stamp && m_particles.size() != 0) {
      lgr::file % "ParticleCounts:" << std::endl;
      for (auto sp : m_particles)
        lgr::file % "\t" << m_properties[sp].name << " = "
                         << m_particles[sp].size() << std::endl;
    }

    // ----- before this line particles are all within borders --- //
    if (stamp) {
      lgr::file % "* ParticleActions" << std::endl;
    }

    {
      auto bundle = typename ParticleAction_t::Bundle_t{
          m_particles, m_J, m_ptc_buffer, m_properties, m_E,  m_B,
          m_grid,      ens, dt,           timestep,     m_rng};
      for (const auto& act : m_ptc_actions) {
        if (!act) continue;
        // TODO check act.range not empty?
        TIMING(
            " -- " + act->name(), START { (*act)(bundle); });
      }
    }

    {  // ----- Before this line _J is local on each cpu, assemble it onto
       // primaries --- //
      TIMING(
          "ReduceJmesh", START {
            // TODOL reduce number of communications?
            for (int i = 0; i < 3; ++i)
              ens.reduce_to_chief(mpi::by::SUM, m_J[i].data().data(),
                                  m_J[i].data().size());
          });
      if (m_cart_opt) {
        TIMING(
            " -- MergeSyncJ",
            START { field::merge_sync_guard_cells(m_J, *m_cart_opt); });
      }
    }

    // export data before m_J is disrupted by field solver
    if (m_sch_export.is_do(timestep)) {
      TIMING(
          "ExportData", START {
            ExportBundle_t bd{m_particles, m_properties, m_E, m_B, m_J,
                              m_grid,      m_cart_opt,   ens, dt,  timestep};
            if (m_f_prior_export) (*m_f_prior_export)(bd);

            for (const auto& exporter : m_exporters) {
              const auto& range = exporter.get_range();
              if (apt::range::is_empty(range)) continue;
              exporter.export_data(bd.timestep, bd.dt, m_sch_export.num_files,
                                   bd.cart_opt, bd.ens, bd.grid, bd.E, bd.B,
                                   bd.J, bd.particles, bd.properties);
            }

            if (m_f_post_export) (*m_f_post_export)(bd);
          });
    }

    if (m_sch_vitals.is_do(timestep)) {
      TIMING(
          "PrintVitals", START {
            print_vitals(m_this_run_dir + "/vitals.txt", timestep * dt);
          });
    }

    if (m_cart_opt) {
      const auto& cart = *m_cart_opt;
      auto bundle = typename FieldAction_t::Bundle_t{
          m_E, m_B, m_J, m_grid, cart, timestep, dt};
      if (stamp) {
        lgr::file % "* FieldActions" << std::endl;
      }

      for (const auto& act : m_fld_actions) {
        // TODO check range empty
        if (!act) continue;
        TIMING(
            " -- " + act->name(), START { (*act)(bundle); });
      }

      TIMING(
          " -- CopySyncEB", START {
            // NOTE sub_range same as domain_range FIXME rethink domain and
            // sub range
            field::copy_sync_guard_cells(m_E, cart);
            field::copy_sync_guard_cells(m_B, cart);
          });
    }
  }

  if (m_sch_dlb.is_do(timestep)) {
    TIMING(
        "DynamicLoadBalance", START {
          // TODO has a few hyper parameters
          // TODO touch create is not multinode safe even buffer is used
          std::optional<int> old_label;
          if (m_ens_opt) old_label.emplace(m_ens_opt->label());

          dynamic_load_balance(m_particles, m_ens_opt, m_cart_opt,
                               m_sch_dlb.target_load);

          std::optional<int> new_label;
          if (m_ens_opt) new_label.emplace(m_ens_opt->label());
          if (stamp) lgr::file % "  update parts of simulator..." << std::endl;
          if (old_label != new_label) update_parts();
          if (stamp) lgr::file % "  Done." << std::endl;
        });
  }

  static ckpt::Autosave autosave;  // significant only on mpi::world.rank() == 0
  if (m_sch_ckpt.is_do(timestep)) {
    TIMING(
        "SaveCheckpoint", START {
          int num_files = m_sch_ckpt.num_files;
          int max_num_checkpoints = m_sch_ckpt.max_num_checkpoints;
          auto dir = ckpt::save_checkpoint(m_this_run_dir, num_files, m_ens_opt,
                                           timestep, m_E, m_B, m_particles,
                                           m_properties);
          if (mpi::world.rank() == 0) {
            auto obsolete_ckpt =
                autosave.add_checkpoint(dir, max_num_checkpoints);
            if (obsolete_ckpt) fs::remove_all(*obsolete_ckpt);
          }
        });
  }

  if (stamp_all) {
    lgr::file.indent_reset();
    lgr::file << "  Total Lapse = " << stamp_all->lapse().in_units_of("ms")
              << std::endl;
  }
}
}  // namespace pic
