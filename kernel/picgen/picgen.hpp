#ifndef  _PICGEN_HPP_
#define  _PICGEN_HPP_

#include "abstract_field_updater.hpp"
#include "kernel/grid.hpp"

#include "old_field_solver/old_field_solver_adapter.hpp"
#include "particle_updater.hpp"
#include "field/communication.hpp"

#include "dye/dynamic_balance.hpp"

#include "particle/migration.hpp"

#include <memory>

#include "gen.hpp"

namespace pic {
  template < int DGrid,
             typename Real,
             template < typename > class PtcSpecs,
             typename ShapeF,
             typename RealJ,
             knl::coordsys CS >
  struct PIC {
  private:
    const knl::Grid< Real, DGrid >& _supergrid;
    const int _guard;
    std::optional<mpi::CartComm> _cart_opt;
    util::Rng<Real> _rng;

    knl::Grid< Real, DGrid > _grid;
    std::optional<dye::Ensemble<DGrid>> _ens_opt;
    apt::array< apt::pair<Real>, DGrid > _borders;

    field::Field<Real, 3, DGrid> _E;
    field::Field<Real, 3, DGrid> _B;
    field::Field<RealJ, 3, DGrid> _J;
    particle::map<particle::array<Real, PtcSpecs>> _particles;

    std::unique_ptr<AbstractFieldUpdater<Real,DGrid,RealJ>> _field_update;
    std::unique_ptr<particle::ParticleUpdater<DGrid,Real,PtcSpecs,ShapeF,RealJ,CS>> _ptc_update;

    std::vector<particle::cParticle<Real, PtcSpecs>> _migrators;

    FieldBC_FoldBackJ<true, DGrid, Real, PtcSpecs, RealJ> _fbj_lower;
    FieldBC_FoldBackJ<false, DGrid, Real, PtcSpecs, RealJ> _fbj_upper;
    Injector< DGrid, Real, PtcSpecs, RealJ> _injector;

    Real _unit_e;

    // ScalarField<Scalar> pairCreationEvents; // record the number of pair creation events in each cell.
    // PairCreationTracker pairCreationTracker;

    constexpr auto all_but( field::offset_t all, int but_comp ) noexcept {
      apt::array<field::offset_t, DGrid> res;
      apt::foreach<0,DGrid>( [&all]( auto& ofs ) { ofs = all; }, res );
      if ( but_comp < DGrid ) res[but_comp] = !all;
      return res;
    }

    void refresh( const dye::Ensemble<DGrid>& ens, double unit_e ) {
      apt::Index<DGrid> bulk_dims;
      for ( int i = 0; i < DGrid; ++i ) {
        int dim = _supergrid[i].dim() / ens.cart_dims[i];
        _grid[i] = _supergrid[i];
        _grid[i].clip( ens.cart_coords[i] * dim, dim );
        // TODOL cart_dim = 1 and periodic
        _borders[i] = { _grid[i].lower(), _grid[i].upper() };
        bulk_dims[i] = _grid[i].dim();
      }

      {
        _E = {{ bulk_dims, _guard }};
        _B = {{ bulk_dims, _guard }};
        // NOTE minimum number of guards of J on one side is ( supp + 1 ) / 2 + 1
        _J = {{ bulk_dims, ( ShapeF::support() + 3 ) / 2 }};

        for( int i = 0; i < 3; ++i ) {
          _E.set_offset( i, all_but( MIDWAY, i ) );
          _B.set_offset( i, all_but( INSITU, i ) );
          _J.set_offset( i, all_but( MIDWAY, i ) );
        }
        _E.reset();
        _B.reset();
        _J.reset();
      }

      ens.is_at_boundary();
      if ( _cart_opt )
        _field_update.reset(new ::ofs::OldFieldUpdater<>( unit_e, *_cart_opt, _grid, ens.is_at_boundary(), _guard ) );
      _ptc_update.reset(new particle::ParticleUpdater<DGrid, Real, PtcSpecs, ShapeF, RealJ, CS>( _grid, _rng ) );
    }

  public:
    PIC( const knl::Grid< Real, DGrid >& supergrid, const std::optional<mpi::CartComm>& cart_opt, int guard, util::Rng<Real> rng, double unit_e )
      : _supergrid(supergrid), _guard(guard), _cart_opt(cart_opt), _rng(std::move(rng)),
        _injector{ _grid, _E, _B, _J, _particles },
        _fbj_lower{ _grid, _E, _B, _J, _particles },
        _fbj_upper{ _grid, _E, _B, _J, _particles },
        _unit_e(unit_e) {
      _grid = supergrid;
      _ens_opt = dye::create_ensemble<DGrid>(cart_opt);
      if ( !_ens_opt ) return;

      const auto& ens = *_ens_opt;
      refresh(ens, unit_e);

      InitialConditionDipole ic( _grid, _E, _B, _J, _particles );
      ic();
    }

    void evolve( int timestep, Real dt ) {
      if ( _ens_opt ) {
        const auto& ens = *_ens_opt;

        // TODOL Opimize communication. Use persistent and buffer?
        for ( int i = 0; i < 3; ++i ) {
          auto& buffer = _J[i].data();
          ens.intra.template reduce<mpi::IN_PLACE>( mpi::by::SUM, ens.chief, buffer.data(), buffer.size() );
        }

        if ( _cart_opt ) {
          field::merge_guard_cells_into_bulk( _J, *_cart_opt );
          (*_field_update)(_E, _B, _J, dt, timestep);
        }

        // TODOL reduce number of communications?
        for ( int i = 0; i < 3; ++i )
          ens.intra.broadcast( ens.chief, _E[i].data().data(), _E[i].data().size() );
        for ( int i = 0; i < 3; ++i )
          ens.intra.broadcast( ens.chief, _B[i].data().data(), _B[i].data().size() );

        _J.reset();
        // if ( false )
        //   sort_particles();
        (*_ptc_update) ( _particles, _J, _E, _B, dt, timestep );

        _fbj_lower();
        _fbj_upper();

        { // migration
          auto migrate_dir =
            []( auto q, auto lb, auto ub ) noexcept {
              return ( q >= lb ) + ( q > ub );
            };

          for ( auto&[ sp, ptcs ] : _particles ) {
            for ( auto ptc : ptcs ) { // TODOL semantics
              char mig_dir = 0;
              for ( int i = 0; i < DGrid; ++i ) {
                mig_dir += migrate_dir( ptc.q()[i], _borders[i][LFT], _borders[i][RGT] ) * apt::pow3(i);
              }

              if ( mig_dir != ( apt::pow3(DGrid) - 1 ) / 2 ) {
                _migrators.emplace_back(std::move(ptc));
                _migrators.back().extra() = mig_dir;
              }
            }
          }

          particle::migrate( _migrators, ens.inter, timestep );
          for ( auto&& ptc : _migrators ) {
            _particles[ptc.template get<particle::species>()].push_back( std::move(ptc) );
          }
          _migrators.resize(0);
        }

        _injector( timestep, dt, _rng );
      }


      // TODO NOTE injection and annihilation are more like free sources and sinks of particles users control, in other words they fit the notion of particle boundary conditions. Do them outside of Updater. In contrast, pair creation is a physical process we model, so it is done inside the updater.
      // inject_particles(); // TODO only coordinate space current is needed in implementing current regulated injection // TODO compatible with ensemble??
      // TODOL annihilation will affect deposition // NOTE one can deposit in the end
      // annihilate_mark_pairs( );

      if (false) {
        // TODO check idle?
        // TODO export_data();
        if ( false ) {
          // TODO has a few hyper parameters
          // TODO touch create is not multinode safe even buffer is used
          std::optional<int> old_label;
          if ( _ens_opt ) old_label.emplace(_ens_opt->label());

          dynamic_load_balance( _particles, _ens_opt, _cart_opt, 100000 );

          std::optional<int> new_label;
          if ( _ens_opt ) new_label.emplace(_ens_opt->label());
          if ( old_label != new_label ) refresh(*_ens_opt, _unit_e);
        }
      }

      // TODOL
      // if (false)
      //   save_snapshot();

      // TODOL
      // if (false)
      //   save_tracing();
    }
  };
}

#endif