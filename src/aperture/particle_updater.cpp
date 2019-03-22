#include "./particle_updater.hpp"
#include "./ensemble.hpp"

#include "particle/pusher.hpp"
#include "particle/pair_producer.hpp"
#include "particle/properties.hpp"
#include "particle/migration.hpp"

#include "field/field_shape_interplay.hpp"
#include <algorithm>

#include "dynamic_variables.hpp"
#include "parameters.hpp"
#include "kernel/coordinate.hpp"

#include "kernel/grid.hpp"
#include "kernel/shapef.hpp"

namespace aperture::impl {
  template < particle::species sp, particle::PairScheme pair_scheme, knl::coordsys CS,
             typename DynaVars_t, typename dJField, typename Migrators, typename Params_t, typename ShapeF, typename Rng, typename Borders >
  void update_species( DynaVars_t& dvars, dJField& dJ, Migrators& migrators, const Params_t& params, const ShapeF& shapef, Rng& rng, const Borders& borders ) {
    using namespace particle;
    if ( !dvars.particles.has(sp) ) return;

    auto dt = params.dt;

    auto charge_over_dt = charge_x<sp> * params.e / params.dt;
    for ( auto& ptc : dvars[sp] ) { // TODOL check ptc is proxy
      if( ptc.is(flag::empty) ) continue;

      if constexpr ( is_charged<sp> ) {
          auto E = field::interpolate(dvars.E, ptc.q, shapef );
          auto B = field::interpolate(dvars.E, ptc.q, shapef );
          auto&& dp = update_p<sp>( ptc, dt, E, B );

          // TODOL make sure the newly added particles will not participate in the loop
          if constexpr ( pair_scheme != PairScheme::Disabled && is_radiative<sp> ) {
            auto Rc = calc_Rc( dt, ptc.p, std::move(dp) );
            auto gamma = std::sqrt( is_massive<sp> + apt::sqabs(ptc.p) );
            if ( !is_productive_lepton( ptc, gamma, Rc, rng ) )
              goto end_of_pair_produce;

            if constexpr ( pair_scheme == PairScheme::Photon ) {
              produce_photons( std::back_inserter( dvars.photons ),
                               ptc, gamma, Rc );
              } else if ( pair_scheme == PairScheme::Instant ) {
              instant_produce_pairs( std::back_inserter( dvars.electrons ),
                                     std::back_inserter( dvars.positrons ),
                                     ptc, gamma, Rc );
            }
          }
          end_of_pair_produce:;

        }
      else if ( sp == species::photon && pair_scheme == PairScheme::Photon ) {
        photon_produce_pairs( std::back_inserter( dvars.electrons ),
                              std::back_inserter( dvars.positrons ),
                              ptc );
      }

      // NOTE q is updated, starting from here, particles may be in the guard cells.
      auto&& dq = update_q<sp,CS>( ptc, dt );
      // TODOL pusher handle boundary condition. Is it needed?
      if constexpr ( is_charged<sp> )
                     dJ.deposit( charge_over_dt, ptc.q()-std::move(dq), ptc.q(), shapef ); // TODOL check the 2nd argument

      if ( is_migrate( ptc.q(), borders ) ) {
        migrators.emplace_back(std::move(ptc));
      }

    }
  }

  template < particle::PairScheme PS, knl::coordsys CS, int I = 0, typename... Args >
  inline void iterate_species( Args&&... args ) {
    if constexpr ( I == 4 ) return; // TODOL loop over all species. Would void_t help?
    else {
      // update_species<static_cast<species>(I), PS, CS>( std::forward<Args>(args)... );
      return iterate_species<PS, CS, I+1>( std::forward<Args>(args)... );
    }
  }

}

namespace aperture {
  template < typename Real, int DGrid, int DPtc, typename state_t,
             knl::shape Shape, particle::PairScheme pair_scheme, knl::coordsys CS
             >
  void ParticleUpdater< Real, DGrid, DPtc, state_t, Shape, pair_scheme, CS >
  ::operator() ( int timestep, DynaVars_t& dvars,
                 const Params_t& params,
                 const std::optional<mpi::CartComm>& cart,
                 const Ensemble<DGrid>& ensemble ) {
    using namespace particle;
    // borders tell which particles will be migrated. These are simply the boundaries of the bulk grid
    apt::array< apt::pair<Real>, DGrid > borders;
    apt::foreach<0,DGrid>
      ( []( auto& b, const auto& g ) {
          b[LFT] = g.lower();
          b[RGT] = g.upper();
        }, borders, dvars.E.mesh().bulk() );

    _dJ.reset();

    // TODOL reduce number of communications?
    for ( int i = 0; i < 3; ++i )
      ensemble.intra.broadcast( ensemble.chief, dvars.E[i].data().data(), dvars.E[i].data().size() );
    for ( int i = 0; i < 3; ++i )
      ensemble.intra.broadcast( ensemble.chief, dvars.B[i].data().data(), dvars.B[i].data().size() );

    impl::iterate_species<pair_scheme, CS>( dvars, _dJ, _migrators, params,
                                            knl::shapef_t<Shape>(), _rng, borders );

    migrate( _migrators, _intercomms, borders, timestep );
    for ( auto&& ptc : _migrators ) {
      dvars.particles[ptc.template get<species>()].push_back( std::move(ptc) ); // TODOL check this
    }
    _migrators.resize(0);

    _dJ.reduce( ensemble.chief, ensemble.intra );

    if ( cart )
      _dJ.integrate( *cart ); // TODO return the deposited J?

  }

}

#include "./traits.hpp"
using namespace traits;
namespace aperture {
  template class ParticleUpdater< real_t, DGrid, DPtc, ptc_state_t, shape,
                                  pair_produce_scheme, coordinate_system>;
}
