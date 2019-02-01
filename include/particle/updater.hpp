#ifndef  _PARTICLE_UPDATER_HPP_
#define  _PARTICLE_UPDATER_HPP_

namespace particle {
  enum class PairScheme
    { Disabled, Instant, Photon };
}

#include "particle/pusher.hpp"
#include "particle/depositer.hpp"
#include "particle/pair_producer.hpp"
#include "particle/migration.hpp"

#include "dynamic_variables.hpp"
#include "parameters.hpp"
#include "field/interpolation.hpp"
#include <algorithm>

namespace todo {
  // NOTE WJ can be in long double, but J only needs to be in double
  void getJ_from_WJ( auto& J, const auto& WJ ) {
    ensemble.reduce(WJ, root);
    if ( is_root ) copy(J, WJ);
    //TODO do the rest of J communication and stuff
  }
}

namespace particle {
  // TODO this is supposed to be used in iterate_species
  constexpr auto present_species;

  template < typename Real_t, std::size_t DGrid, std::size_t DPtc,
             sf::shape S, PairScheme pair_scheme, CoordSys CS,
             species posion // posion = positron || ion in injection
             >
  class Updater {
  private:
    using WJ_t = Real_t;
    field::Field< WJ_t, 3, DGrid > _WJ;
    vector<Real_t, DPtc> _migrators;
    Rng _rng;

    template < species sp, typename Dynvar_t, typename Params_t, typename Grid >
    inline void update_species( Dynvar_t& dvars, const Params_t& params, const Grid& grid ) {
      auto dt = params.dt;

      for ( auto& ptc : dvars.particles<sp> ) {
        if( ptc.is(flag::empty) ) continue;

        if constexpr ( is_charged<sp> ) {
          auto E = field::interpolate(dvars.E, grid, ptc.q);
          auto B = field::interpolate(dvars.E, grid, ptc.q);
          auto&& dp = update_p<sp>( ptc, dt, E, B );

          // TODO make sure the newly added particles will not participate in the loop
          if constexpr ( pair_scheme != PairScheme::Disabled && is_radiative<sp> ) {
            auto Rc = calc_Rc( dt, ptc.p, std::move(dp) );
            auto gamma = std::sqrt( is_massive<sp> + apt::abs_sq(ptc.p) );
            if ( !is_productive_lepton( ptc, gamma, Rc, _rng ) )
              goto end_of_pair_produce;

            if constexpr ( pair_scheme == PairScheme::Photon ) {
              produce_photons( std::back_inserter( dvars.particles<species::photon> ),
                               ptc, gamma, Rc );
            } else if ( pair_scheme == PairScheme::Instant ) {
              instant_produce_pairs( std::back_inserter( dvars.particles<species::electron> ),
                                     std::back_inserter( dvars.particles<species::positron> ),
                                     ptc, gamma, Rc );
            }
          }
          end_of_pair_produce:

        } else if ( sp == species::photon && pair_scheme == PairScheme::Photon ) {
          photon_produce_pairs( std::back_inserter( dvars.particles<species::electron> ),
                                std::back_inserter( dvars.particles<species::positron> ),
                                ptc );
        }

        // NOTE q is updated, starting from here, particles may be in the guard cells.
        auto&& dq = update_q<sp,CS>( ptc, dt );
        // pusher handle boundary condition
        if constexpr ( is_charged<sp> ) depositWJ<S>( _WJ, ptc, std::move(dq), grid );

        if ( is_migrate( ptc.q, params.borders ) ) {
          // TODO make sure after move, the moved from ptc is set to flag::empty
          _migrators.push_back(std::move(ptc));
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

    template < std::size_t I = 0, typename... Args >
    inline void iterate_species( Args&... args ) {
      if constexpr ( I == std::tuple_size_v<decltype(present_species)> ) return;
      update_species<std::get<I>(present_species)>( args... );
      return iterate_species<I+1>( args... );
    }


  public:
    void operator() ( DynamicVars<Real_t, DGrid, DPtc>& dvars,
                      const Params<Real_t, DGrid>& params,
                      const Grid<DGrid>& grid, auto& comm ) {
      const auto& borders; // TODO
      const auto& neighbors = params.neighbors;

      for ( auto& comp : _WJ ) for ( auto& elm : comp ) elm = 0.0;
      EnsBroadcastFields( dvars.E, dvars.B ); // TODO

      //-------------------------------
      inject( fetch( species::electron, dvars ), fetch( posion, dvars ) );

      // NOTE one can deposit in the end
      annihilate_mark_pairs( );

      //-------------------------------
      // TODO use compile time known
      iterate_species( dvars, params, grid );

      migrate( _migrators, params.neighbors, borders, comm );
      for ( auto&& ptc : _migrators ) {
        fetch(ptc.get<species>(), dvars).push_back( std::move(ptc) );
      }
      _migrators.resize(0);

      // do this before J is actually used
      todo::getJ_from_WJ( dvars.J, _WJ );
    }

  } update;

}

#endif
