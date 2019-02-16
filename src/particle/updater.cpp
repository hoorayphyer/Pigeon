#include "particle/updater.hpp"

#include "particle/pusher.hpp"
#include "particle/depositer.hpp"
#include "particle/pair_producer.hpp"
#include "particle/properties.hpp"

#include "field/interpolation.hpp"
#include <algorithm>

#include "dynamic_variables.hpp"
#include "parameters.hpp"
#include "kernel/coordinate.hpp"

#include "kernel/grid.hpp"
#include "kernel/shapef.hpp"

namespace particle {
  template < typename Real, int DGrid, int DPtc, typename state_t, knl::shape Shape,
             PairScheme pair_scheme, knl::coordsys CS, species posion >
  template < species sp, typename Dynvar_t, typename Params_t, typename Grid >
  void Updater<Real,DGrid,DPtc,state_t,Shape,pair_scheme,CS,posion>
  ::update_species( Dynvar_t& dvars, const Params_t& params, const Grid& grid ) {
      if ( dvars[sp].size() == 0 ) return;

      auto dt = params.dt;

      for ( auto&& ptc : dvars[sp] ) {
        if( ptc.is(flag::empty) ) continue;

        if constexpr ( is_charged<sp> ) {
          auto E = field::interpolate(dvars.E, grid, ptc.q);
          auto B = field::interpolate(dvars.E, grid, ptc.q);
          auto&& dp = update_p<sp>( ptc, dt, E, B );

          // TODO make sure the newly added particles will not participate in the loop
          if constexpr ( pair_scheme != PairScheme::Disabled && is_radiative<sp> ) {
            auto Rc = calc_Rc( dt, ptc.p, std::move(dp) );
            auto gamma = std::sqrt( is_massive<sp> + apt::sqabs(ptc.p) );
            if ( !is_productive_lepton( ptc, gamma, Rc, _rng ) )
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

        } else if ( sp == species::photon && pair_scheme == PairScheme::Photon ) {
          photon_produce_pairs( std::back_inserter( dvars.electrons ),
                                std::back_inserter( dvars.positrons ),
                                ptc );
        }

        // NOTE q is updated, starting from here, particles may be in the guard cells.
        auto&& dq = update_q<sp,CS>( ptc, dt );
        // pusher handle boundary condition
        if constexpr ( is_charged<sp> ) depositWJ<Shape>( _WJ, ptc, std::move(dq), grid );

        if ( is_migrate( ptc.q, params.borders ) ) {
          // TODO make sure after move, the moved from ptc is set to flag::empty
          _migrators.emplace_back(std::move(ptc));
          // ptc.set(flag::empty);
        }

      }

      if constexpr ( is_charged<sp> ) {
        auto tmp = -1 * charge_x<sp> * params.e * grid[0].delta / dt;
        for ( auto& elm : _WJ[0] ) elm *= tmp;

        tmp *= (grid[1].delta / grid[0].delta);
        for ( auto& elm : _WJ[1] ) elm *= tmp;

        if constexpr ( DGrid == 2 ) {
          tmp = charge_x<sp> * params.e;
        } else if ( DGrid == 3 ) {
          tmp *= (grid[2].delta / grid[1].delta);
        } static_assert(DGrid < 4);

        for ( auto& elm : _WJ[2] ) elm *= tmp;

        // // TODO check here if needs to output species J
        // if ( is_dataexport ) {
        //   data_out.WJ[sp] = WJ;
        //   todo::getJ_from_WJ(data_out.J[sp], data_out.WJ[sp] );
        // }
      }

      // TODO sort particle array periodically
    }

  template < typename Real, int DGrid, int DPtc, typename state_t, knl::shape Shape,
             PairScheme pair_scheme, knl::coordsys CS, species posion >
  void Updater<Real,DGrid,DPtc,state_t,Shape,pair_scheme,CS,posion>
  ::operator() ( DynamicVars<Real, DGrid, DPtc, state_t>& dvars,
                 const Params<Real, DGrid>& params,
                 const knl::Grid<DGrid, Real>& grid,
                 const mpi::Comm& ensemble ) {
      const auto& borders = params.borders;
      const auto& neighbors = params.neighbors;

      for ( auto& comp : _WJ ) for ( auto& elm : comp ) elm = 0.0;
      // EnsBroadcastFields( dvars.E, dvars.B ); // TODO
      //-------------------------------
      inject( dvars.electrons, dvars[posion] );

      // NOTE one can deposit in the end
      // annihilate_mark_pairs( ); // TODO

      //-------------------------------
      // TODO use compile time known
      iterate_species( dvars, params, grid );

      migrate( _migrators, params.neighbors, borders, _intercomms );
      for ( auto&& ptc : _migrators ) {
        // TODO something wrong
        // dvars[ptc.get<species>()].push_back( std::move(ptc) );
      }
      _migrators.resize(0);

      // // TODO do this before J is actually used
      // todo::getJ_from_WJ( dvars.J, _WJ );
    }

}

#include "traits.hpp"
using namespace traits;
namespace particle {
  template class Updater< real_t, DGrid, DPtc, ptc_state_t, shape,
                          pair_produce_scheme,
                          coordinate_system, posion_inj >;
}
