#ifndef  _APERTURE_HPP_
#define  _APERTURE_HPP_

#include "abstract_field_updater.hpp"
#include "abstract_particle_updater.hpp"
#include "kernel/grid.hpp"
#include "kernel/shapef.hpp"

#include "old_field_solver/old_field_solver_adapter.hpp"
#include "particle_updater.hpp"

#include "ensemble/dynamic_balance.hpp"

#include "traits.hpp"
#include <memory>

namespace aperture::tmp {
  using namespace traits;
  template < typename Real, int DGrid, typename state_t >
  using PtcUpdater = ParticleUpdater<Real, DGrid, 3, state_t, knl::shapef_t<shape>, pair_produce_scheme, coordinate_system >;
}

namespace aperture {
  template < typename Real, int DGrid, typename state_t >
  struct Aperture {
  private:
    std::optional<mpi::CartComm> _cart_opt;
    std::optional<Ensemble<DGrid>> _ens_opt;
    knl::Grid< Real, DGrid, knl::grid1d::Clip > _grid;

    field::Field<Real, 3, DGrid> _E;
    field::Field<Real, 3, DGrid> _B;
    field::Field<Real, 3, DGrid> _J;
    particle::map<particle::array<Real,3,state_t>> _particles;

    std::unique_ptr<AbstractFieldUpdater<Real,DGrid>> _field_update;
    std::unique_ptr<AbstractParticleUpdater<Real, DGrid, state_t>> _ptc_update;

    // ScalarField<Scalar> pairCreationEvents; // record the number of pair creation events in each cell.
    // PairCreationTracker pairCreationTracker;

  public:
    Aperture( const std::optional<mpi::CartComm>& cart_opt, const knl::Grid< Real, DGrid >& supergrid, int guard )
      : _cart_opt(cart_opt), _grid(knl::make_clip(supergrid)), _E({_grid, guard}), _B({_grid, guard}), _J({_grid, guard}) {
      // TODO check local grid constructor. Set anchor and dim

      // create cart and ensemble, set local grid anchor
      _ens_opt = create_ensemble<DGrid>(cart_opt);
      if ( !_ens_opt ) return;
      const auto& ens = *_ens_opt;
      if ( _cart_opt )
        _field_update.reset(new ofs::OldFieldUpdater<>( *_cart_opt, _grid, ens.is_at_boundary(), guard ) );
      util::Rng<Real> rng; // TODO set seed
      _ptc_update.reset(new tmp::PtcUpdater<Real, DGrid, state_t>( _grid, rng, _cart_opt, ens ) );
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
        (*_ptc_update)( _J, _particles, _E, _B, dt, unit_e, timestep );
      }

      if (false) {
        // TODO check idle?
        // TODO export_data();
        if ( false ) {
          // TODO has a few hyper parameters
          // TODO also update localgrid, and maybe fieldupdater and ptcupdater as well
          // TODO touch create is not multinode safe even buffer is used
          dynamic_load_balance( _particles, _ens_opt, _cart_opt, 100000 );
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
