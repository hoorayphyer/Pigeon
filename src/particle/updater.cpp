#include "particle/updater.hpp"

#include "particle/pusher.hpp"
#include "particle/pair_producer.hpp"
#include "particle/properties.hpp"

#include "field/field_shape_interplay.hpp"
#include <algorithm>

#include "dynamic_variables.hpp"
#include "parameters.hpp"
#include "kernel/coordinate.hpp"

#include "kernel/grid.hpp"
#include "kernel/shapef.hpp"

namespace particle::impl {
  template < species sp, PairScheme pair_scheme, knl::coordsys CS,
             typename DynaVars_t, typename WJField, typename Migrators, typename Params_t, typename ShapeF, typename Rng, typename Borders >
  void update_species( DynaVars_t& dvars, WJField& WJ, Migrators& migrators, const Params_t& params, const ShapeF& shapef, Rng& rng, const Borders& borders ) {
    if ( dvars[sp].size() == 0 ) return;

    auto dt = params.dt;

    for ( auto& ptc : dvars[sp] ) {
      if( ptc.is(flag::empty) ) continue;

      if constexpr ( is_charged<sp> ) {
          auto E = field::interpolate(dvars.E, ptc.q, shapef );
          auto B = field::interpolate(dvars.E, ptc.q, shapef );
          auto&& dp = update_p<sp>( ptc, dt, E, B );

          // TODO make sure the newly added particles will not participate in the loop
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
      // pusher handle boundary condition
      if constexpr ( is_charged<sp> )
        field::depositWJ( WJ, charge_x<sp>, ptc.q(), std::move(dq), shapef );

      if ( is_migrate( ptc.q(), borders ) ) {
        migrators.emplace_back(std::move(ptc));
      }

    }

    // TODO sort particle array periodically
  }

  template < PairScheme PS, knl::coordsys CS, int I = 0, typename... Args >
  inline void iterate_species( Args&&... args ) {
    if constexpr ( I == 4 ) return; // TODOL loop over all species. Would void_t help?
    else {
      // update_species<static_cast<species>(I), PS, CS>( std::forward<Args>(args)... );
      return iterate_species<PS, CS, I+1>( std::forward<Args>(args)... );
    }
  }

}

namespace particle {
  template < typename Real, int DGrid, int DPtc, typename state_t,
             knl::shape Shape, PairScheme pair_scheme, knl::coordsys CS
             >
  void Updater< Real, DGrid, DPtc, state_t, Shape, pair_scheme, CS >
  ::operator() ( int timestep, DynaVars_t& dvars,
                 const Params_t& params,
                 const std::optional<mpi::Comm>& ensemble ) {

    // TODO borders at real physical boundary need to be specified still, because now mesh controls margin
    apt::array< apt::pair<Real>, DGrid > borders; // TODO borders not done yet

    apt::foreach<0, 3>
      ( []( auto comp ) { // returns a proxy
          for ( auto& elm : comp.data() ) elm = 0.0;
        }, _WJ );

      // EnsBroadcastFields( dvars.E, dvars.B ); // TODO
      //-------------------------------
      // inject( dvars.electrons, dvars[posion] ); // TODO only coordinate space current is needed in implementing current regulated injection

      // NOTE one can deposit in the end
      // annihilate_mark_pairs( ); // TODO

      //-------------------------------
    // TODO NOTE one can use one chunk of memory for WJ so that only one pass of reduce is needed
    // TODO NOTE injection and annihilation are more like free sources and sinks of particles users control, in other words they fit the notion of particle boundary conditions. Do them outside of Updater. In contrast, pair creation is a physical process we model, so it is done inside the updater.
    impl::iterate_species<pair_scheme, CS>( dvars, _WJ, _migrators, params,
                             knl::shapef_t<Shape>(), _rng, borders );

    migrate( _migrators, _intercomms, borders, timestep );
    for ( auto&& ptc : _migrators ) {
      dvars.particles[ptc.template get<species>()].push_back( std::move(ptc) ); // TODO check this
    }
    _migrators.resize(0);

    // TODO do this before J is actually used
    {
      const auto& grid = dvars.E.mesh().bulk();

      apt::foreach<0, DGrid>
        ( [&]( auto comp, const auto& g ) {
            auto tmp = -1 * params.e * g.delta() / params.dt;
            for ( auto& elm : comp.data() ) elm *= tmp;
          }, _WJ, grid );

      if constexpr ( DGrid == 2 ) {
        auto tmp = params.e / params.dt;
        for ( auto& elm : _WJ[2].data() ) elm *= tmp;
      }

      // Maybe one can use
      // TODO getJ_from_WJ( dvars.J, _WJ );
      // }
      // TODOL now WJ mixes all species, how to output species components of the current?
    }
  }

}

#include "traits.hpp"
using namespace traits;
namespace particle {
  template class Updater< real_t, DGrid, DPtc, ptc_state_t, shape,
                          pair_produce_scheme, coordinate_system>;
}
