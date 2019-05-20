#ifndef  _PIC_SIMULATOR_HPP_
#define  _PIC_SIMULATOR_HPP_

#include "manifold/grid.hpp"

#include "field/updater.hpp"
#include "particle/updater.hpp"
#include "field/communication.hpp"

#include "dye/dynamic_balance.hpp"

#include "particle/migration.hpp"

#include "io/io.hpp"

#include <memory>

#include "gen.hpp"

#include "logger/logger.hpp"

namespace pic {
  std::string this_run_dir;


  template < int DGrid,
             typename Real,
             template < typename > class PtcSpecs,
             typename ShapeF,
             typename RealJ,
             typename Metric >
  struct Simulator {
  private:
    const mani::Grid< Real, DGrid >& _supergrid;
    const int _guard;
    std::optional<mpi::CartComm> _cart_opt;
    util::Rng<Real> _rng;

    mani::Grid< Real, DGrid > _grid;
    std::optional<dye::Ensemble<DGrid>> _ens_opt;
    apt::array< apt::pair<Real>, DGrid > _borders;

    field::Field<Real, 3, DGrid> _E;
    field::Field<Real, 3, DGrid> _B;
    field::Field<RealJ, 3, DGrid> _J;
    particle::map<particle::array<Real, PtcSpecs>> _particles;

    std::unique_ptr<field::Updater<Real,DGrid,RealJ>> _field_update;
    std::unique_ptr<particle::Updater<DGrid,Real,PtcSpecs,ShapeF,RealJ,Metric>> _ptc_update;

    std::vector<particle::cParticle<Real, PtcSpecs>> _migrators;

    std::unique_ptr<bc::Axissymmetric<DGrid, Real, PtcSpecs, RealJ>> _fbc_axis;
    std::unique_ptr<bc::FoldBackJ<DGrid, Real, PtcSpecs, RealJ>> _fbj;
    std::unique_ptr<Injector< DGrid, Real, PtcSpecs, RealJ>> _injector;

    // ScalarField<Scalar> pairCreationEvents; // record the number of pair creation events in each cell.
    // PairCreationTracker pairCreationTracker;

    constexpr auto all_but( field::offset_t all, int but_comp ) noexcept {
      apt::array<field::offset_t, DGrid> res;
      apt::foreach<0,DGrid>( [&all]( auto& ofs ) { ofs = all; }, res );
      if ( but_comp < DGrid ) res[but_comp] = !all;
      return res;
    }

    void refresh( const dye::Ensemble<DGrid>& ens ) {
      apt::Index<DGrid> bulk_dims;
      for ( int i = 0; i < DGrid; ++i ) {
        _grid[i] = _supergrid[i].divide( ens.cart_dims[i], ens.cart_coords[i] );
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

      _fbc_axis.reset( new typename decltype(_fbc_axis)::element_type {_grid, _E, _B, _J, _particles} );
      _fbj.reset( new typename decltype(_fbj)::element_type {_grid, _E, _B, _J, _particles} );
      _injector.reset( new typename decltype(_injector)::element_type {_grid, _E, _B, _J, _particles} );

      ens.is_at_boundary();
      if ( _cart_opt )
        _field_update.reset(new field::Updater<Real,DGrid,RealJ>( *_cart_opt, _grid, ens.is_at_boundary(), _guard ) );
      _ptc_update.reset(new particle::Updater<DGrid, Real, PtcSpecs, ShapeF, RealJ, Metric>( _grid, _rng ) );
    }

  public:
    Simulator( const mani::Grid< Real, DGrid >& supergrid, const std::optional<mpi::CartComm>& cart_opt, int guard )
      : _supergrid(supergrid), _guard(guard), _cart_opt(cart_opt) {
      _grid = supergrid;
      _ens_opt = dye::create_ensemble<DGrid>(cart_opt);
      if ( !_ens_opt ) return;

      const auto& ens = *_ens_opt;
      refresh(ens);
    }

    template < template < int, typename, template < typename > class, typename > class IC >
    int load_initial_condition() {
      // TODOL a temporary fix, which may crash under the edge case in which initially many particles are created
      int init_ts = 0;
      if (_cart_opt) {
        IC ic( _grid, _E, _B, _J, _particles );
        ic();
        init_ts = ic.initial_timestep();
      }
      // broadcast to ensure uniformity
      mpi::world.broadcast( 0, &init_ts, 1 );
      return init_ts;
    }

    inline void set_rng_seed( int seed ) { _rng.set_seed(seed); }

    void evolve( int timestep, Real dt ) {
      if ( _ens_opt ) {
        const auto& ens = *_ens_opt;

        // lgr::file % "reduce J" << std::endl;
        // TODOL Opimize communication. Use persistent and buffer?
        for ( int i = 0; i < 3; ++i ) {
          auto& buffer = _J[i].data();
          ens.intra.template reduce<mpi::IN_PLACE>( mpi::by::SUM, ens.chief, buffer.data(), buffer.size() );
        }

        if ( _cart_opt ) {
          // lgr::file % "merge guard cells of J" << std::endl;
          field::merge_guard_cells_into_bulk( _J, *_cart_opt );
          // lgr::file % "field update" << std::endl;
          (*_field_update)(_E, _B, _J, dt, timestep);
          // lgr::file % "field BC" << std::endl;
          (*_fbc_axis)();
        }

        // lgr::file % "broadcast E and B" << std::endl;
        // TODOL reduce number of communications?
        for ( int i = 0; i < 3; ++i )
          ens.intra.broadcast( ens.chief, _E[i].data().data(), _E[i].data().size() );
        for ( int i = 0; i < 3; ++i )
          ens.intra.broadcast( ens.chief, _B[i].data().data(), _B[i].data().size() );

        _J.reset();
        // if ( false )
        //   sort_particles();

        // lgr::file % "particle update" << std::endl;
        (*_ptc_update) ( _particles, _J, _E, _B, dt, timestep );

        // lgr::file % "particle BC" << std::endl;
        (*_fbj)();

        // lgr::file % "migration" << std::endl;
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

        // lgr::file % "injection" << std::endl;
        // (*_injector)( timestep, dt, _rng );
      }


      // TODO NOTE injection and annihilation are more like free sources and sinks of particles users control, in other words they fit the notion of particle boundary conditions. Do them outside of Updater. In contrast, pair creation is a physical process we model, so it is done inside the updater.
      // inject_particles(); // TODO only coordinate space current is needed in implementing current regulated injection // TODO compatible with ensemble??
      // TODOL annihilation will affect deposition // NOTE one can deposit in the end
      // annihilate_mark_pairs( );

      if ( timestep >= pic::data_export_init_ts && (timestep % pic::interval::data_export == 0 ) && _ens_opt ) {
        // lgr::file % "export_data" << std::endl;
        io::export_data<pic::real_export_t, pic::DGrid, pic::real_t, particle::Specs, pic::ShapeF, pic::real_j_t, pic::Metric>( this_run_dir, timestep, dt, pic::pmpio_num_files, _cart_opt, *_ens_opt, _grid, _E, _B, _J, _particles  );
        if ( false ) {
          // TODO has a few hyper parameters
          // TODO touch create is not multinode safe even buffer is used
          std::optional<int> old_label;
          if ( _ens_opt ) old_label.emplace(_ens_opt->label());

          // TODO manually initialize all species in the simulation
          dynamic_load_balance( _particles, _ens_opt, _cart_opt, 100000 );

          std::optional<int> new_label;
          if ( _ens_opt ) new_label.emplace(_ens_opt->label());
          if ( old_label != new_label ) refresh(*_ens_opt);
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
