#ifndef  _APERTURE_HPP_
#define  _APERTURE_HPP_

#include "abstract_field_updater.hpp"
#include "kernel/grid.hpp"

#include "old_field_solver/old_field_solver_adapter.hpp"
#include "particle_updater.hpp"
#include "field/communication.hpp"

#include "dye/dynamic_balance.hpp"

#include "particle/migration.hpp"

#include <memory>

#include "traits.hpp"
namespace aperture {
  using namespace traits;
  template < typename Real, int DGrid, typename state_t >
  using ParticleUpdater = particle::ParticleUpdater< Real, DGrid, 3, state_t, ShapeF, real_dj_t, coordinate_system >;
}

namespace aperture {
  template < typename Real, int DGrid, typename state_t >
  struct Aperture {
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
    field::Field<Real, 3, DGrid> _J;
    particle::map<particle::array<Real,3,state_t>> _particles;

    std::unique_ptr<AbstractFieldUpdater<Real,DGrid>> _field_update;
    std::unique_ptr<ParticleUpdater<Real, DGrid, state_t>> _ptc_update;

    std::vector<particle::cParticle<Real, DPtc, state_t>> _migrators;

    // ScalarField<Scalar> pairCreationEvents; // record the number of pair creation events in each cell.
    // PairCreationTracker pairCreationTracker;

    void refresh( const dye::Ensemble<DGrid>& ens ) {
      field::Mesh<2> mesh ( knl::dims(_grid), _guard );
      for ( int i = 0; i < DGrid; ++i ) {
        int dim = _supergrid[i].dim() / ens.cart_dims[i];
        _grid[i] = _supergrid[i];
        _grid[i].clip( ens.cart_coords[i] * dim, dim );
        _E = { mesh };
        _B = { mesh };
        _J = { mesh };
        // TODOL cart_dim = 1 and periodic
        _borders[i] = { _grid[i].lower(), _grid[i].upper() };
      }

      if ( _cart_opt )
        _field_update.reset(new ofs::OldFieldUpdater<>( *_cart_opt, _grid, ens.is_at_boundary(), _guard ) );
      _ptc_update.reset(new ParticleUpdater<Real, DGrid, state_t>( _grid, _rng ) );
    }

  public:
    Aperture( const knl::Grid< Real, DGrid >& supergrid, const std::optional<mpi::CartComm>& cart_opt, int guard, util::Rng<Real> rng ) : _supergrid(supergrid), _guard(guard), _cart_opt(cart_opt), _rng(std::move(rng))
    {
      _grid = supergrid;
      _ens_opt = dye::create_ensemble<DGrid>(cart_opt);
      if ( !_ens_opt ) return;

      const auto& ens = *_ens_opt;
      refresh(ens);
    }

    void evolve( int timestep, Real dt, Real unit_e ) {
      if ( _ens_opt ) {
        const auto& ens = *_ens_opt;

        if ( _cart_opt ) (*_field_update)(_E, _B, _J, dt, timestep);
        // TODOL reduce number of communications?
        for ( int i = 0; i < 3; ++i )
          ens.intra.broadcast( ens.chief, _E[i].data().data(), _E[i].data().size() );
        for ( int i = 0; i < 3; ++i )
          ens.intra.broadcast( ens.chief, _B[i].data().data(), _B[i].data().size() );

        // if ( false )
        //   sort_particles();
        auto& Jmesh = (*_ptc_update)( _particles, _E, _B, dt, unit_e, timestep );

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

        // TODO Injection here

        // TODOL Opimize communication. Use persistent and buffer? NOTE although probably not important, reduce might still cause loss of precision, so we perform reduce on Jmesh instead of first assignning Jmesh to _J then reducing on _J
        for ( int i = 0; i < 3; ++i ) {
          auto& buffer = Jmesh[i].data();
          ens.intra.template reduce<mpi::IN_PLACE>( mpi::by::SUM, ens.chief, buffer.data(), buffer.size() );
        }

        if ( _cart_opt ) {
          field::merge_guard_cells_into_bulk( Jmesh, *_cart_opt );
          // Copy Jmesh into J NOTE this extra copy is because Jmesh might be in long double where _J is in double

          const auto mesh_size = Jmesh.mesh().stride(DGrid);
          for ( int iJ = 0; iJ < 3; ++iJ ) {
            auto& dataJ = _J[iJ].data();
            const auto& dataJmesh = Jmesh[iJ].data();
            for ( int i = 0; i < mesh_size; ++i ) dataJ[i] = dataJmesh[i];
          }
        }
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
