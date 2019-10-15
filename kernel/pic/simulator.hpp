#ifndef  _PIC_SIMULATOR_HPP_
#define  _PIC_SIMULATOR_HPP_

#include "manifold/grid.hpp"

#include "particle/updater.hpp"
#include "field/sync.hpp"
#include "field/yee.hpp"

#include "particle/migration.hpp"
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
#include "pic/action_reserve.hpp"

#ifdef PIC_DEBUG
#include "debug/debugger.hpp"
#endif

namespace pic {
  std::string this_run_dir;

  template < int DGrid,
             typename R,
             template < typename > class S,
             typename ShapeF,
             typename RJ >
  struct Simulator {
  private:
    const mani::Grid< R, DGrid >& _supergrid;
    std::optional<mpi::CartComm> _cart_opt;
    util::Rng<R> _rng;

    mani::Grid< R, DGrid > _grid;
    std::optional<dye::Ensemble<DGrid>> _ens_opt;
    apt::array< apt::pair<R>, DGrid > _borders;

    field::Field<R, 3, DGrid> _E;
    field::Field<R, 3, DGrid> _B;
    field::Field<RJ, 3, DGrid> _J;
    particle::map<particle::array<R, S>> _particles;

    std::vector<std::unique_ptr<field::Action<R,DGrid,RJ>>> _field_actions;
    std::vector<std::unique_ptr<particle::Action<DGrid,R,S,RJ>>> _ptc_actions;

    ActionReserve<DGrid,R> _rsv;

    std::vector<particle::Particle<R, S>> _migrators;

    template < typename Action >
    void taylor( Action& a ) {
      // NOTE range is assumed to be noempty [,)
      // POLEDANCE check the logic here. including one more cell at upper boundary
      auto f =
        [] ( int Ib_global, int Ie_global,
             const mani::Grid1D<R>& supergrid,
             const mani::Grid1D<R>& localgrid,
             bool is_periodic,
             int guard
             ) noexcept -> apt::pair<int> {
         // lb is the Ib_global with respect to the localgrid
         int lb = Ib_global - static_cast<int>( ( localgrid.lower() - supergrid.lower() ) / localgrid.delta() + 0.5 );
         int ub = Ie_global - static_cast<int>( ( localgrid.lower() - supergrid.lower() ) / localgrid.delta() + 0.5 );

         // Extend to include guard cells on true boundaries. Whether or not it will be used is controled by Ib_global and Ie_global
         int lb_local = 0;
         int ub_local = localgrid.dim();
         if ( !is_periodic ) {
           if ( std::abs( localgrid.lower() - supergrid.lower() ) < supergrid.delta() )
             lb_local -= guard;
           if ( std::abs( localgrid.upper() - supergrid.upper() ) < supergrid.delta() )
             ub_local += guard;
         }

         if ( ub <= lb_local || lb >= ub_local ) {
           // range not applicable on current local patch
           return {0,0};
         } else {
           lb = std::max<int>( lb_local, lb );
           ub = std::min<int>( ub, ub_local );
           return { lb, ub };
         }
        };

      for ( int i = 0; i < DGrid; ++i ) {
        auto [b_new, e_new] = f( a[i].begin(), a[i].end(), _supergrid[i], _grid[i], pic::periodic[i], field::myguard );
        a[i].begin() = b_new;
        a[i].end() = e_new;
      }
    }

    apt::pair<std::optional<field::Field<R, 3, DGrid>>> backupEB( int i_this_action ) const {
      const auto& fa = _field_actions;
      if ( !fa[i_this_action] ) return {};
      if ( i_this_action >= fa.size() - 1 ) return {};
      if ( !fa[i_this_action+1] || !fa[i_this_action+1]->require_original_EB ) return {};

      const auto& ath = *fa[i_this_action];
      const auto& ane = *fa[i_this_action + 1];

      apt::Index<DGrid> ovIb, ovExt;
      bool is_empty = false;
      for ( int i = 0; i < DGrid; ++i ) {
        ovIb[i] = std::max(ath.actIb[i] - ath.actGuard[i][LFT], ane.actIb[i] - ane.actGuard[i][LFT]);
        ovExt[i] = std::min(ath.actIb[i] + ath.actExt[i] + ath.actGuard[i][RGT], ane.actIb[i] + ane.actExt[i] + ane.actGuard[i][RGT]) - ovIb[i];
        if (ovExt[i] <= 0) is_empty = true;
      }

      if ( is_empty ) return {};

      field::Field<R, 3, DGrid> E_bak({ovExt,0}), B_bak({ovExt,0});
      for ( const auto& I : apt::Block(ovExt) ) E_bak[0](I) = _E[0](I+ovIb);
      for ( const auto& I : apt::Block(ovExt) ) E_bak[1](I) = _E[1](I+ovIb);
      for ( const auto& I : apt::Block(ovExt) ) E_bak[2](I) = _E[2](I+ovIb);

      for ( const auto& I : apt::Block(ovExt) ) B_bak[0](I) = _B[0](I+ovIb);
      for ( const auto& I : apt::Block(ovExt) ) B_bak[1](I) = _B[1](I+ovIb);
      for ( const auto& I : apt::Block(ovExt) ) B_bak[2](I) = _B[2](I+ovIb);

      return {{std::move(E_bak)}, {std::move(B_bak)}};
    }

    void update_parts( const dye::Ensemble<DGrid>& ens ) {
      for ( int i = 0; i < DGrid; ++i ) {
        _grid[i] = _supergrid[i].divide( ens.cart_dims[i], ens.cart_coords[i] );
        // TODO cart_dim = 1 and periodic
        _borders[i] = { _grid[i].lower(), _grid[i].upper() };
      }

      { // set field actions and reserve
        const auto f_actions = field::set_up_field_actions<DGrid,R,RJ>();
        _field_actions.resize(f_actions.size());
        for ( int i = 0; i < f_actions.size(); ++i ) {
          _field_actions[i].reset(f_actions[i]->Clone());
          taylor(*_field_actions[i]);
        }
        // TODO reorder by require_old_EB

        { // only reserve for necessary actions. For example, when there is only one action, or the action doesn't require old_EB, skip it.
          _rsv = {};
          std::vector<const apt::ActionBase<DGrid>*> ptrs;
          int count = 0;
          for ( int i = 0; i < _field_actions.size(); ++i ) {
            if ( !_field_actions[i] || !(_field_actions[i] -> orig_EB()) ) continue;
            ptrs.emplace_back(_field_actions[i].get());
            ++count;
          }
          if ( count > 1) _rsv.init(ptrs);
        }
      }

      { // set particle actions
        const auto p_actions = particle::set_up_particle_actions<DGrid,R,S,ShapeF,RJ>();
        _ptc_actions.resize(p_actions.size());
        for ( int i = 0; i < p_actions.size(); ++i ) {
          _ptc_actions[i].reset(p_actions[i]->Clone());
          taylor(*_ptc_actions[i]);
        }
      }

      if ( pic::msperf_qualified(_ens_opt) ) lgr::file.open(std::ios_base::app);
    }

    void migrate_particles( int timestep ) {
      using namespace particle;
      // bulk range = [lb, ub)
      constexpr auto migrate_code =
        []( auto q, auto lb, auto ub ) noexcept {
          return ( q >= lb ) + ( q >= ub );
        };

      for ( auto sp : _particles ) {
        for ( auto ptc : _particles[sp] ) { // TODOL semantics
          if ( !ptc.is(flag::exist) ) continue;
          migrInt<DGrid> mig_dir{};
          for ( int i = 0; i < DGrid; ++i ) {
            mig_dir += migrate_code( ptc.q()[i], _borders[i][LFT], _borders[i][RGT] ) * apt::pow3(i);
          }

          if ( mig_dir != ( apt::pow3(DGrid) - 1 ) / 2 ) {
            mig_dir.imprint(ptc);
            _migrators.emplace_back(std::move(ptc));
          }
        }
      }

      migrate( _migrators, _ens_opt->inter, timestep );

      for ( auto&& ptc : _migrators ) {
        if ( !ptc.is(flag::exist) ) continue;
        auto sp = ptc.template get<species>();
#ifdef PIC_DEBUG
        // // check if the received ptc trully resides in this ensemble.
        // apt::array<int,DGrid> mig_co;
        // bool is_OK = true;
        // for ( int i = 0; i < DGrid; ++i ) {
        //   mig_co[i] = migrate_code( ptc.q()[i], _borders[i][LFT], _borders[i][RGT] );
        //   if ( mig_co[i] != 1 && !_ens_opt->is_at_boundary(i)[(mig_co[i] != 0)] ) //NOTE need to consider boundaries
        //     is_OK = false;
        // }
        // if ( !is_OK ) {
        //   lgr::file << "ts=" << debug::timestep << ", wr=" << debug::world_rank << ", el=" << debug::ens_label << std::endl;
        //   lgr::file << "Received across-ensemble particles! q = " << ptc.q() << ", p = " << ptc.p() << ", birth = " << static_cast<int>(ptc.template get<birthplace>()) << std::endl;
        //   lgr::file << "  mig_dir on new ensemble  = " << mig_co;
        //   // get old mig_co
        //   for ( int i = 0; i < DGrid; ++i ) {
        //     mig_co[i] = ( migrInt<DGrid>(ptc) % apt::pow3(i+1) ) / apt::pow3(i);
        //   }
        //   lgr::file << ", mig_dir on old ensemble = " << mig_co << std::endl;
        //   debug::throw_error("Received across-ensemble particles!");
        // }
#endif
        ptc.template reset<destination>();
        _particles[sp].push_back( std::move(ptc) );
      }
      _migrators.resize(0);
    }

    template < typename MR >
    inline bool is_do( const MR& mr, int timestep ) const noexcept {
      return mr.is_on && timestep >= mr.init_ts && (timestep % mr.interval == 0 );
    }

  public:
    Simulator( const mani::Grid< R, DGrid >& supergrid, const std::optional<mpi::CartComm>& cart_opt )
      : _supergrid(supergrid), _cart_opt(cart_opt) {
      _grid = supergrid;
      if ( pic::dlb_init_replica_deploy ) {
        std::optional<int> label{};
        if ( _cart_opt )
          label.emplace( _cart_opt->rank() );
        else {
          int num_ens = 1;
          for ( int i = 0; i < DGrid; ++i ) num_ens *= pic::dims[i];

          int my_rank = mpi::world.rank();
          int count = num_ens;
          for ( int i = 0; i < num_ens; ++i ) {
            int num_replicas = (*pic::dlb_init_replica_deploy)(i);
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
        auto range = apt::make_range({}, bulk_dims, field::myguard);
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
      for ( auto sp : particle::properties )
        _particles.insert( sp, particle::array<R, S>() );

      if ( _ens_opt ) update_parts(*_ens_opt);
    }

    int load_initial_condition( std::optional<std::string> checkpoint_dir ) {
      int init_ts = 0;
      if ( checkpoint_dir ) {
        init_ts = ckpt::load_checkpoint( *checkpoint_dir, _ens_opt, _cart_opt, _E, _B, _particles, particle::N_scat );
        ++init_ts; // checkpoint is saved at the end of a timestep
        if ( _ens_opt ) update_parts(*_ens_opt);
      } else {
        // TODOL a temporary fix, which may crash under the edge case in which initially many particles are created
        if (_cart_opt) {
          auto ic = pic::set_up_initial_conditions<DGrid,R,RJ,S>();
          taylor(ic);
          ic(_grid, _E, _B, _J, _particles);
          init_ts = ic.initial_timestep();
          field::copy_sync_guard_cells(_E, *_cart_opt);
          field::copy_sync_guard_cells(_B, *_cart_opt);
          field::copy_sync_guard_cells(_J, *_cart_opt);
        }
      }
      // broadcast to ensure uniformity
      mpi::world.broadcast( 0, &init_ts, 1 );
      // set up vitals
      if ( mpi::world.rank() == 0 ) {
        vital::t_phys_prev = ( init_ts - 1 ) * pic::dt; // NOTE the -1
        for ( auto sp : _particles ) {
          double num = 0.0;
          for ( const auto& ptc : _particles[sp] ) {
            if ( ptc.is(particle::flag::exist) ) num += ptc.frac();
          }
          vital::num_ptcs_prev.push_back(num);
          vital::num_scat_prev.push_back(particle::N_scat[sp]);
        }
      }
      return init_ts;
    }

    inline void set_rng_seed( int seed ) { _rng.set_seed(seed); }

    // NOTE we choose to do particle update before field update so that in saving snapshots we don't need to save _J
    void evolve( int timestep, R dt ) {
      if ( timestep % pic::cout_ts_interval == 0 && mpi::world.rank() == 0 )
        std::cout << "==== Timestep " << timestep << " ====" << std::endl;

#ifdef PIC_DEBUG
      debug::timestep = timestep;
      debug::world_rank = mpi::world.rank();
      if ( _ens_opt ) debug::ens_label = _ens_opt -> label();
#endif

      std::optional<tmr::Timestamp> stamp_all;
      std::optional<tmr::Timestamp> stamp;
      if ( is_do(pic::msperf_mr, timestep) && pic::msperf_qualified(_ens_opt) ) {
        stamp_all.emplace();
        stamp.emplace();
        if ( msperf_max_entries &&
             ( timestep - msperf_mr.init_ts > *msperf_max_entries * msperf_mr.interval )  ) lgr::file.clear();
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

        if ( is_do(pic::sort_particles_mr, timestep) ) {
          TIMING("SortParticles", START {
              for ( auto sp : _particles ) particle::sort( _particles[sp] );
            });
        }

        if ( stamp && _particles.size() != 0 ) {
          lgr::file % "ParticleCounts:" << std::endl;
          for ( auto sp : _particles )
            lgr::file % "\t" << particle::properties[sp].name << " = " << _particles[sp].size() << std::endl;
        }

        // ----- before this line particles are all within borders --- //
        TIMING("ParticleActions", START {
            for ( int i = 0; i < _ptc_actions.size(); ++i ) {
              if ( !_ptc_actions[i] ) continue;
              if ( stamp ) {
                lgr::file % "--" << _ptc_actions[i]->name() << std::endl;
              }
              auto& act = *_ptc_actions[i];
              act( _particles, _J, particle::properties, _E, _B, _grid, &ens, dt, timestep, _rng );
            }
          });

        TIMING("MigrateParticles", START {
            migrate_particles( timestep );
          });

        // ----- After this line, particles are all within borders.

        // ----- Before this line _J is local on each cpu --- //
        TIMING("ReduceJmesh", START {
            // TODOL reduce number of communications?
            for ( int i = 0; i < 3; ++i )
              ens.reduce_to_chief( mpi::by::SUM, _J[i].data().data(), _J[i].data().size() );
          });

        if ( _cart_opt ) {
          TIMING("FieldUpdate", START {
              field::merge_sync_guard_cells( _J, *_cart_opt );
              _rsv.reserve(_E,_B);
              for ( int i = 0; i < _field_actions.size(); ++i ) {
                if ( !_field_actions[i] ) continue;
                if ( stamp ) {
                  lgr::file % "--" << _field_actions[i]->name() << std::endl;
                }
                const auto& act = *_field_actions[i];
                // FIXME
                // std::cout << "ts = " << timestep << std::endl;
                // std::cout << _field_actions[i]->name() << std::endl;
                // std::cout << "  revert" << std::endl;
                // _rsv.revert_to_prior(_E,_B,act);
                // std::cout << "  action" << std::endl;
                act(_E, _B, _J, _grid, *_cart_opt, timestep, dt);
                // std::cout << "  back" << std::endl;
                // _rsv.back_to_current(_E,_B);
              }
            });
        }
      }

      if ( is_do(pic::export_data_mr, timestep) && _ens_opt ) {
        TIMING("ExportData", START {
            auto fexps = io::set_up_field_export<real_export_t,DGrid,real_t,real_j_t>();
            auto pexps = io::set_up_particle_export<real_export_t,DGrid,real_t,particle::Specs>();

            io::set_is_collinear_mesh(io::is_collinear_mesh); // TODO

            io::export_data<pic::real_export_t>( this_run_dir, timestep, dt, pic::pmpio_num_files, pic::downsample_ratio, _cart_opt, *_ens_opt, _grid, _E, _B, _J, _particles, fexps, pexps );
            for ( auto ptr : fexps ) delete ptr;
            for ( auto ptr : pexps ) delete ptr;
          });
      }

      if ( is_do(pic::dlb_mr, timestep) ) {
        TIMING("DynamicLoadBalance", START {
            if ( _ens_opt ) { // first reduce N_scat to avoid data loss
              const auto& ens = *_ens_opt;
              ens.reduce_to_chief( mpi::by::SUM, particle::N_scat.data().data(), particle::N_scat.data().size() );
              if ( !ens.is_chief() ) {
                for ( auto& x : particle::N_scat.data() ) x = 0;
              }
            }
            // TODO has a few hyper parameters
            // TODO touch create is not multinode safe even buffer is used
            std::optional<int> old_label;
            if ( _ens_opt ) old_label.emplace(_ens_opt->label());

            dynamic_load_balance( _particles, _ens_opt, _cart_opt, pic::dlb_target_load );

            std::optional<int> new_label;
            if ( _ens_opt ) new_label.emplace(_ens_opt->label());
            if (stamp) lgr::file % "  update parts of simulator..." << std::endl;
            if ( old_label != new_label ) update_parts(*_ens_opt);
            if (stamp) lgr::file % "  Done." << std::endl;
          });
      }

      if (_ens_opt && is_do(pic::vitals_mr, timestep) ) {
        TIMING("Statistics", START {
            pic::check_vitals( pic::this_run_dir + "/vitals.txt", timestep * dt, *_ens_opt, _cart_opt, _particles, particle::N_scat );
          });
      }

      static ckpt::Autosave autosave; // significant only on mpi::world.rank() == 0
      if ( is_do(pic::checkpoint_mr, timestep)
           || ( pic::checkpoint_autosave_hourly &&
                autosave.is_save({*pic::checkpoint_autosave_hourly * 3600, "s"}) ) ) {
        TIMING("SaveCheckpoint", START {
            if ( _ens_opt ) { // first reduce N_scat to avoid data loss
              const auto& ens = *_ens_opt;
              ens.reduce_to_chief( mpi::by::SUM, particle::N_scat.data().data(), particle::N_scat.data().size() );
              if ( !ens.is_chief() ) {
                for ( auto& x : particle::N_scat.data() ) x = 0;
              }
            }
            auto dir = ckpt::save_checkpoint( this_run_dir, num_checkpoint_parts, _ens_opt, timestep, _E, _B, _particles, particle::N_scat );
            if ( mpi::world.rank() == 0 ) {
              auto obsolete_ckpt = autosave.add_checkpoint(dir, pic::max_num_ckpts);
              if ( obsolete_ckpt ) fs::remove_all(*obsolete_ckpt);
            }
          });
      }

      // TODOL
      // if (false)
      //   save_tracing();
      if ( stamp_all ) {
        lgr::file.indent_reset();
        lgr::file << "  Total Lapse = " << stamp_all->lapse().in_units_of("ms") << std::endl;
      }
    }
  };
}

#endif
