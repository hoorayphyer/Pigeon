#include "particle/particle_updater.hpp"
#include "particle/pusher.hpp"
#include "particle/depositer.hpp"
#include "particle/pair_producer.hpp"
#include "particle/migration.hpp"

namespace particle {
  namespace todo {
    // NOTE WJ can be in long double, but J only needs to be in double
    void getJ_from_WJ( auto& J, const auto& WJ ) {
      ensemble.reduce(WJ, root);
      if ( is_root ) copy(J, WJ);
      //TODO do the rest of J communication and stuff
    }
  }

  template < sf::shape S, PairScheme pair_scheme, CoordSys CS >
  void update( ) {
    auto dt;
    auto& rng;
    auto& WJ;
    auto& comm;
    auto& migrate_buffer; // for all particles
    const auto& bounds;
    const auto& neighbors;

    WJ.zero_out();
    EnsBroadcastFields( Efield, Bfield );

    auto& data;
    auto& electrons, positrons;
    //-------------------------------
    // posion = positron || ion
    inject( electrons, posion );

    // NOTE one can deposit in the end
    annihilate_mark_pairs( );

    //-------------------------------
    // TODO use compile time known
    static_for ( auto[sp, particles] : particles_map ) {

      for ( auto& ptc : particles ) {
        if( ptc.is(flag::empty) ) continue;

        if constexpr ( is_charged<sp> ) {
          auto E, B; // TODO interpolate field
          auto&& dp = update_p<sp>( ptc, sp, dt, E, B );

          // TODO make sure the newly added particles will not participate in the loop
          if constexpr ( pair_scheme != PairScheme::Disabled && is_radiative<sp> ) {
            auto Rc = calc_Rc( dt, ptc.p, std::move(dp) );
            auto gamma = std::sqrt( (sp.mass > 0) + apt::abs_sq(ptc.p) );
            if ( !is_productive_lepton( ptc, gamma, Rc, rng ) )
              goto end_of_pair_produce;

            if constexpr ( pair_scheme == PairScheme::Photon ) {
              produce_photons( std::back_inserter( data[species::photon] ),
                               ptc, gamma, Rc );
            } else if ( pair_scheme == PairScheme::Instant ) {
              instant_produce_pairs( std::back_inserter( data[species::electron] ),
                                     std::back_inserter( data[species::positron] ),
                                     ptc, gamma, Rc );
            }
          }
          end_of_pair_produce:

        } else if ( sp == species::photon && pair_scheme == PairScheme::Photon ) {
          photon_produce_pairs( std::back_inserter( data[species::electron] ),
                                std::back_inserter( data[species::positron] ),
                                ptc );
        }

        // NOTE q is updated, starting from here, particles may be in the guard cells.
        auto&& dq = update_q<sp,CS>( ptc, sp, dt );
        // pusher handle boundary condition
        if constexpr ( is_charged<sp> ) depositWJ<S>( WJ, ptc, std::move(dq), grid );

        if ( is_migrate( ptc.q, bounds ) ) {
          // TODO make sure after move, the moved from ptc is set to flag::empty
          migrate_buffer.push_back(std::move(ptc));
          // ptc.set(flag::empty);
        }

      }


      if constexpr ( is_charged<sp> ) {
        WJ *= charge<sp> * -grid.delta / dt;

        // check here if needs to output species J
        if ( is_dataexport ) {
          data_out.WJ[sp] = WJ;
          todo::getJ_from_WJ(data_out.J[sp], data_out.WJ[sp] );
        }
      }

      // TODO sort particle array periodically

    }

    migrate( migrate_buffer, neighbors, bounds, comm );
    for ( auto&& ptc : migrate_buffer ) {
      particles[ptc.get<species>()].push_back( std::move(ptc) );
    }
    migrate_buffer.resize();

    // do this before J is actually used
    auto& J;
    todo::getJ_from_WJ(J, WJ);

  }
}
