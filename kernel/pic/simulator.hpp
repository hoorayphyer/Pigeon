#ifndef  _PIC_SIMULATOR_HPP_
#define  _PIC_SIMULATOR_HPP_

#include "particle/updater.hpp"
#include "field/sync.hpp"
#include "field/yee.hpp"

#include "particle/sorter.hpp"
#include "io/io.hpp"
#include "dye/dynamic_balance.hpp"
#include "ckpt/checkpoint.hpp"
#include "ckpt/autosave.hpp"

#include <memory>

#include "pic_impl.hpp"

#include "logger/logger.hpp"
#include "timer/timer.hpp"
#include "filesys/filesys.hpp"

#include "pic/vitals.hpp"

#ifdef PIC_DEBUG
#include "debug/debugger.hpp"
#endif

#include <unordered_map>

namespace pic {
  std::string this_run_dir;

  template < int DGrid,
             typename R,
             template < typename > class S,
             typename ShapeF,
             typename RJ >
  struct Simulator {
  private:
    const apt::Grid< R, DGrid >& _supergrid;
    std::optional<mpi::CartComm> _cart_opt;
    util::Rng<R> _rng;

    apt::Grid< R, DGrid > _grid;
    std::optional<dye::Ensemble<DGrid>> _ens_opt;

    particle::map<particle::Properties> _properties;
    field::Field<R, 3, DGrid> _E;
    field::Field<R, 3, DGrid> _B;
    field::Field<RJ, 3, DGrid> _J;
    particle::map<particle::array<R, S>> _particles;

    std::vector<std::unique_ptr<field::Action<R,DGrid,RJ>>> _field_actions;
    std::vector<std::unique_ptr<particle::Action<DGrid,R,S,RJ>>> _ptc_actions;

    std::vector<particle::Particle<R, S>> _ptc_buffer;

    void taylor( apt::ActionBase<DGrid>& a ) {
      // NOTE range is assumed to be noempty [,)
      auto f =
        [] ( int Ib_global, int Ie_global,
             const apt::Grid1D<R>& supergrid,
             const apt::Grid1D<R>& localgrid,
             bool is_periodic
             ) noexcept -> apt::pair<int>
        {
         // shift is the index of local lower in the global grid
         int shift = static_cast<int>( ( localgrid.lower() - supergrid.lower() ) / localgrid.delta() + 0.5 );
         Ib_global -= shift;
         Ie_global -= shift;
         // now Ib_ and Ie_global are with respect to the current local grid

         // Extend to infinity on true boundaries
         // NOTE there may be actions done completely in the guard cells, such as assigining values
         int lb_local = 0;
         int ub_local = localgrid.dim();
         if ( !is_periodic ) {
           if ( std::abs( localgrid.lower() - supergrid.lower() ) < supergrid.delta() )
             lb_local = std::numeric_limits<int>::min();
           if ( std::abs( localgrid.upper() - supergrid.upper() ) < supergrid.delta() )
             ub_local = std::numeric_limits<int>::max();
         }

         if ( Ie_global <= lb_local || Ib_global >= ub_local ) {
           // range not applicable on current local patch
           return {0,0};
         } else {
           Ib_global = std::max<int>( lb_local, Ib_global );
           Ie_global = std::min<int>( Ie_global, ub_local );
           return { Ib_global, Ie_global };
         }
        };

      for ( int i = 0; i < DGrid; ++i ) {
        auto [b_new, e_new] = f( a[i].begin(), a[i].end(), _supergrid[i], _grid[i], pic::periodic[i] );
        a[i].begin() = b_new;
        a[i].end() = e_new;
      }
    }

    void update_parts( const dye::Ensemble<DGrid>& ens ) {
      for ( int i = 0; i < DGrid; ++i ) {
        _grid[i] = _supergrid[i].divide( ens.cart_topos[i].dim(), ens.cart_coords[i] );
      }

      { // set field actions
        // auto f_actions = field::set_up_field_actions<DGrid,R,RJ>();
        auto f_actions = set_up_field_actions();
        _field_actions.resize(f_actions.size());
        for ( int i = 0; i < f_actions.size(); ++i ) {
          _field_actions[i].reset(f_actions[i]->Clone());
          taylor(*_field_actions[i]);
        }
      }

      { // set particle actions
        // auto p_actions = particle::set_up_particle_actions<DGrid,R,S,ShapeF,RJ>();
        auto p_actions = set_up_particle_actions();
        _ptc_actions.resize(p_actions.size());
        for ( int i = 0; i < p_actions.size(); ++i ) {
          _ptc_actions[i].reset(p_actions[i]->Clone());
          taylor(*_ptc_actions[i]);
        }
      }

      { // init runtime data
        RTD::data().init( _properties, _grid );
      }

      if ( profiling_plan.on and profiling_plan.is_qualified() ) lgr::file.open(std::ios_base::app);
    }

  public:
    Simulator( const apt::Grid< R, DGrid >& supergrid, const std::optional<mpi::CartComm>& cart_opt, const particle::map<particle::Properties>& props )
      : _supergrid(supergrid), _cart_opt(cart_opt), _properties(props) {
      _grid = supergrid;
      if ( pic::init_replica_deploy ) {
        std::optional<int> label{};
        if ( _cart_opt )
          label.emplace( _cart_opt->rank() );
        else {
          int num_ens = 1;
          for ( int i = 0; i < DGrid; ++i ) num_ens *= pic::dims[i];

          int my_rank = mpi::world.rank();
          int count = num_ens;
          for ( int i = 0; i < num_ens; ++i ) {
            int num_replicas = (*pic::init_replica_deploy)(i);
            if ( my_rank >= count && my_rank < count + num_replicas ) {
              label.emplace(i);
              break;
            } else
              count += num_replicas;
          }
        }

        auto intra_opt = mpi::world.split(label);
        _ens_opt = dye::create_ensemble<DGrid>(cart_opt, intra_opt);
      } else {
        _ens_opt = dye::create_ensemble<DGrid>(cart_opt);
      }

      { // initialize fields. Idles also do this for the sake of loading checkpoint, so pic::dims is used instead of _ens_opt or _cart_opt
        apt::Index<DGrid> bulk_dims;
        for ( int i = 0; i < DGrid; ++i ) {
          bulk_dims[i] = _supergrid[i].dim() / pic::dims[i];
        }
        // FIXME : should include upper boundary ???
        auto range = apt::make_range({}, bulk_dims, myguard);
        _E = {range};
        _B = {range};
        _J = {range};

        for( int i = 0; i < 3; ++i ) {
          _E.set_offset( i, field::yee::ofs_gen<DGrid>( field::yee::Etype, i ) );
          _B.set_offset( i, field::yee::ofs_gen<DGrid>( field::yee::Btype, i ) );
          _J.set_offset( i, field::yee::ofs_gen<DGrid>( field::yee::Etype, i ) );
        }
        _E.reset();
        _B.reset();
        _J.reset();
      }
      // NOTE all species in the game should be created regardless of whether they appear on certain processes. This is to make the following work
      // 1. detailed balance. Absence of some species may lead to deadlock to transferring particles of that species.
      // 2. data export. PutMultivar requires every patch to output every species
      for ( auto sp : _properties )
        _particles.insert( sp, particle::array<R, S>() );

      if ( _ens_opt ) update_parts(*_ens_opt);
    }

    int load_initial_condition( std::optional<std::string> checkpoint_dir ) {
      int init_ts = 0;
      if ( checkpoint_dir ) {
        init_ts = ckpt::load_checkpoint( *checkpoint_dir, _ens_opt, _cart_opt, _E, _B, _particles, _properties );
        ++init_ts; // checkpoint is saved at the end of a timestep
        if ( _ens_opt ) update_parts(*_ens_opt);
      } else {
        // FIXME a temporary fix, which may crash under the edge case in which initially many particles are created.
        if (_cart_opt) {
          auto ic = set_up_initial_conditions();
          taylor(ic);
          ic(_grid, _E, _B, _J, _particles);
          for ( auto sp : _particles ) {
            for ( auto ptc : _particles[sp] ) { // TODOL semantics
              ptc.set(sp);
            }
          }
          field::copy_sync_guard_cells(_E, *_cart_opt);
          field::copy_sync_guard_cells(_B, *_cart_opt);
          field::copy_sync_guard_cells(_J, *_cart_opt);
        }
      }
      // broadcast to ensure uniformity
      mpi::world.broadcast( 0, &init_ts, 1 );

      // initialize trace_counters to ensure unique trace serial numbers across runs
      ParticleTracing::init(_properties, _particles);

      return init_ts;
    }

    inline void set_rng_seed( int seed ) { _rng.set_seed(seed); }

    // NOTE we choose to do particle update before field update so that in saving snapshots we don't need to save _J
    void evolve( int timestep, R dt ) {
      if ( timestep % print_timestep_to_stdout_interval == 0 && mpi::world.rank() == 0 )
        std::cout << "==== Timestep " << timestep << " ====" << std::endl;

#ifdef PIC_DEBUG
      debug::timestep = timestep;
      debug::world_rank = mpi::world.rank();
      if ( _ens_opt ) debug::ens_label = _ens_opt -> label();
#endif

      std::optional<tmr::Timestamp> stamp_all;
      std::optional<tmr::Timestamp> stamp;
      if ( profiling_plan.is_do(timestep) and profiling_plan.is_qualified() ) {
        stamp_all.emplace();
        stamp.emplace();
        if ( profiling_plan.is_reached_max_entries(timestep) ) lgr::file.clear();
        lgr::file << "==== Timestep " << timestep << " ====" << std::endl;
        lgr::file.indent_append("\t");
      }

      auto TIMING = [&stamp](std::string message, auto f) {
                      if ( stamp ) {
                        lgr::file % message << "==>>" << std::endl;
                        stamp.emplace();
                      }
                      f();
                      if ( stamp ) {
                        lgr::file % "\tLapse " << stamp->lapse().in_units_of("ms") << std::endl;
                      }
                    };
#define START [&]()

      if ( _ens_opt ) {
        const auto& ens = *_ens_opt;

        TIMING("BroadcastEB", START {
            // TODOL reduce number of communications?
            for ( int i = 0; i < 3; ++i )
              ens.intra.broadcast( ens.chief, _E[i].data().data(), _E[i].data().size() );
            for ( int i = 0; i < 3; ++i )
              ens.intra.broadcast( ens.chief, _B[i].data().data(), _B[i].data().size() );
          });

        _J.reset();

        if ( sort_ptcs_plan.is_do(timestep) ) {
          TIMING("SortParticles", START {
              for ( auto sp : _particles ) particle::sort( _particles[sp], _grid );
            });
        }

        if ( stamp && _particles.size() != 0 ) {
          lgr::file % "ParticleCounts:" << std::endl;
          for ( auto sp : _particles )
            lgr::file % "\t" << _properties[sp].name << " = " << _particles[sp].size() << std::endl;
        }

        // ----- before this line particles are all within borders --- //
        if ( stamp ) {
          lgr::file % "* ParticleActions" << std::endl;
        }
        for ( int i = 0; i < _ptc_actions.size(); ++i ) {
          if ( !_ptc_actions[i] ) continue;
          TIMING(" -- " + _ptc_actions[i]->name(), START {
            auto& act = *_ptc_actions[i];
            act( _particles, _J, &_ptc_buffer, _properties, _E, _B, _grid, &ens, dt, timestep, _rng );
          });
        }

        { // NOTE rescale Jmesh back to real grid delta
          auto dV = apt::dV(_grid);

          for ( int i = 0; i < DGrid; ++i ) {
            R tmp = _grid[i].delta() / dV;
            for ( auto& elm : _J[i].data() ) elm *= tmp;
          }

          for ( int i = DGrid; i < 3; ++i ) {
            for ( auto& elm : _J[i].data() ) elm /= dV;
          }
        }

        // ----- Before this line _J is local on each cpu --- //
        TIMING("ReduceJmesh", START {
            // TODOL reduce number of communications?
            for ( int i = 0; i < 3; ++i )
              ens.reduce_to_chief( mpi::by::SUM, _J[i].data().data(), _J[i].data().size() );
          });

        if ( _cart_opt ) {
          if ( stamp ) {
            lgr::file % "* FieldActions" << std::endl;
          }
          TIMING(" -- MergeSyncJ", START {
              field::merge_sync_guard_cells( _J, *_cart_opt );
            });

          for ( int i = 0; i < _field_actions.size(); ++i ) {
            if ( !_field_actions[i] ) continue;
            TIMING(" -- " + _field_actions[i]->name(), START {
              const auto& act = *_field_actions[i];
              act(_E, _B, _J, _grid, *_cart_opt, timestep, dt);
            });
          }

          TIMING(" -- CopySyncEB", START {
              // NOTE sub_range same as domain_range FIXME rethink domain and sub range
              field::copy_sync_guard_cells(_E, *_cart_opt );
              field::copy_sync_guard_cells(_B, *_cart_opt );
            });
        }
      }

      if ( export_plan.is_do(timestep) and _ens_opt ) {
        TIMING("ExportData", START {
            export_prior_hook( _particles, _properties, _E, _B, _J, _grid, *_ens_opt, dt, timestep );

            auto fexps = set_up_field_export();
            auto pexps = set_up_particle_export();

            io::set_is_collinear_mesh(is_collinear_mesh); // TODO

            io::export_data<pic::real_export_t>( this_run_dir, timestep, dt,
                                                 export_plan.num_files,
                                                 export_plan.downsample_ratio,
                                                 _cart_opt, *_ens_opt, _grid, _E, _B, _J,
                                                 _particles, _properties, fexps, pexps );
            for ( auto ptr : fexps ) delete ptr;
            for ( auto ptr : pexps ) delete ptr;

            export_post_hook();
          });
      }

      if ( load_balance_plan.is_do(timestep) ) {
        TIMING("DynamicLoadBalance", START {
            // TODO has a few hyper parameters
            // TODO touch create is not multinode safe even buffer is used
            std::optional<int> old_label;
            if ( _ens_opt ) old_label.emplace(_ens_opt->label());

            dynamic_load_balance( _particles, _ens_opt, _cart_opt, load_balance_plan.target_load );

            std::optional<int> new_label;
            if ( _ens_opt ) new_label.emplace(_ens_opt->label());
            if (stamp) lgr::file % "  update parts of simulator..." << std::endl;
            if ( old_label != new_label ) update_parts(*_ens_opt);
            if (stamp) lgr::file % "  Done." << std::endl;
          });
      }

      if (_ens_opt and vitals_plan.is_do(timestep) ) {
        TIMING("Statistics", START {
            pic::check_vitals( pic::this_run_dir + "/vitals.txt", timestep * dt, *_ens_opt, _cart_opt, _properties, _particles, RTD::data().N_scat );
          });
      }

      static ckpt::Autosave autosave; // significant only on mpi::world.rank() == 0
      if ( checkpoint_plan.is_do(timestep) ) {
        TIMING("SaveCheckpoint", START {
            auto dir = ckpt::save_checkpoint( this_run_dir, checkpoint_plan.num_files, _ens_opt, timestep, _E, _B, _particles, _properties );
            if ( mpi::world.rank() == 0 ) {
              auto obsolete_ckpt = autosave.add_checkpoint(dir, checkpoint_plan.max_num_checkpoints);
              if ( obsolete_ckpt ) fs::remove_all(*obsolete_ckpt);
            }
          });
      }

      if ( save_tracing_plan.is_do(timestep) ) { // FIXME
        TIMING("SaveTracing", START {
            auto dir = ckpt::save_tracing( this_run_dir, save_tracing_plan.num_files, _ens_opt, timestep, _particles );
          });
      }

      if ( stamp_all ) {
        lgr::file.indent_reset();
        lgr::file << "  Total Lapse = " << stamp_all->lapse().in_units_of("ms") << std::endl;
      }
    }
  };
}

#endif
