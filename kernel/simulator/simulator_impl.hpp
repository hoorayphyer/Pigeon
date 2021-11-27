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
  // TODO also taylor exporters

  // TODO
  // {  // init runtime data
  //   RTD::data().init(m_properties, m_grid);
  // }
}

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
void Simulator<DGrid, R, S, RJ, RD>::initialize(
    apt::Grid<R, DGrid>&& supergrid, std::optional<mpi::CartComm>&& cart_opt,
    particle::map<particle::Properties>&& props,
    apt::array<bool, DGrid>&& periodic, const apt::array<int, DGrid>& dims) {
  m_supergrid = std::move(supergrid);
  m_cart_opt = std::move(cart_opt);
  m_properties = std::move(props);
  m_periodic = std::move(periodic);

  m_grid = m_supergrid;

  if (m_sch_prof.on and m_sch_prof.is_qualified()) {
    lgr::file.open(std::ios_base::app);
    lgr::file << "--- Initializing Simulator..." << std::endl;
  }

  // TODO maybe move this to builder. init_replica_deploy
  if (false) {  // if (pic::init_replica_deploy) {
    // std::optional<int> label{};
    // if (m_cart_opt)
    //   label.emplace(m_cart_opt->rank());
    // else {
    //   int num_ens = 1;
    //   for (int i = 0; i < DGrid; ++i) num_ens *= dims[i];

    //   int my_rank = mpi::world.rank();
    //   int count = num_ens;
    //   for (int i = 0; i < num_ens; ++i) {
    //     int num_replicas = (*pic::init_replica_deploy)(i);
    //     if (my_rank >= count && my_rank < count + num_replicas) {
    //       label.emplace(i);
    //       break;
    //     } else
    //       count += num_replicas;
    //   }
    // }

    // auto intra_opt = mpi::world.split(label);
    // m_ens_opt = dye::create_ensemble<DGrid>(m_cart_opt, intra_opt);
  } else {
    m_ens_opt = dye::create_ensemble<DGrid>(m_cart_opt);
  }

  {  // initialize fields. Idles also do this for the sake of loading
     // checkpoint, so dims is used instead of m_ens_opt or m_cart_opt
    apt::Index<DGrid> bulk_dims;
    for (int i = 0; i < DGrid; ++i) {
      bulk_dims[i] = m_supergrid[i].dim() / dims[i];
    }
    auto range = apt::make_range({}, bulk_dims, m_fld_guard);
    m_E = {range};
    m_B = {range};
    m_J = {range};

    for (int i = 0; i < 3; ++i) {
      m_E.set_offset(i, field::yee::ofs_gen<DGrid>(field::yee::Etype, i));
      m_B.set_offset(i, field::yee::ofs_gen<DGrid>(field::yee::Btype, i));
      m_J.set_offset(i, field::yee::ofs_gen<DGrid>(field::yee::Etype, i));
    }
    m_E.reset();
    m_B.reset();
    m_J.reset();
  }
  // NOTE all species in the game should be created regardless of whether they
  // appear on certain processes. This is to make the following work
  // 1. detailed balance. Absence of some species may lead to deadlock to
  // transferring particles of that species.
  // 2. data export. PutMultivar requires every patch to output every species
  for (auto sp : m_properties) m_particles.insert(sp, particle::array<R, S>());

  if (m_ens_opt) update_parts();

  if (m_sch_prof.on and m_sch_prof.is_qualified()) {
    lgr::file << "--- Done." << std::endl;
  }
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
            ExportBundle_t bundle{m_particles, m_properties, m_E,        m_B,
                                  m_J,         m_grid,       m_cart_opt, ens,
                                  dt,          timestep};
            if (m_f_prior_export) (*m_f_prior_export)(bundle);

            // TODO sort out exporter vs action
            // for (auto& exporter : set_up_data_exporters()) {
            //   auto range = exporter.get_range();
            //   taylor(range);
            //   if (apt::range::is_empty(range)) continue;
            //   exporter.export_data(timestep, dt, export_plan.num_files,
            //                        _cart_opt, ens, _grid, _E, _B, _J,
            //                        _particles, _properties);
            // }

            if (m_f_post_export) (*m_f_post_export)(bundle);
          });
    }

    if (m_f_custom_step) {
      // TODO what should its bundle be
      (*m_f_custom_step)();
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
