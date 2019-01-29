#include "particle_module/particle_updater.hpp"
#inlcude "particle_module/pusher.hpp"
#include "particle_module/depositer.hpp"
#include "particle_module/pair_producer.hpp"
#include "particle_module/migration.hpp"

namespace particle {
  namespace todo {
    auto find_Rc( auto dp ) {
    }

    // NOTE WJ can be in long double, but J only needs to be in double
    void getJ_from_WJ( auto& J, const auto& WJ ) {
      ensemble.reduce(WJ, root);
      if ( is_root ) copy(J, WJ);
      // do the rest of J communication and stuff
    }
  }

  template < sf::shape S, pair_scheme PS, CoordSys CS >
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
        if( ptc.is<flag::empty>() ) continue;

        if constexpr ( sp.is_charged ) {
          auto E, B; // TODO interpolate field
          auto&& dp = update_p( ptc, sp, dt, E, B );

          // TODO make sure the newly added particles will not participate in the loop
          if constexpr ( PS != DISABLED && sp.is_radiative ) {
            auto Rc = calc_Rc( dt, ptc.p, std::move(dp) );
            auto gamma = std::sqrt( (sp.mass > 0) + apt::abs_sq(ptc.p) );
            if ( !is_productive_lepton( ptc, gamma, Rc, rng ) )
              goto end_of_pair_produce;

            if constexpr ( PS == PHOTON ) {
              produce_photons( std::back_inserter( std::get<Photon>(data) ),
                               ptc, gamma, Rc );
            } else if ( PS == INSTANT ) {
              instant_produce_pairs( std::back_inserter( std::get<Electron>(data) ),
                                     std::back_inserter( std::get<Positron>(data) ),
                                     ptc, gamma, Rc );
            }
          }
          end_of_pair_produce:

        } else if ( sp == photon && PS == PHOTON ) {
          photon_produce_pairs( std::back_inserter( std::get<Electron>(data) ),
                                std::back_inserter( std::get<Positron>(data) ),
                                ptc );
        }

        // NOTE q is updated, starting from here, particles may be in the guard cells.
        auto&& dq = update_q<CS>( ptc, sp, dt );
        // pusher handle boundary condition
        if constexpr ( sp.is_charged ) depositWJ<S>( WJ, ptc, std::move(dq), grid );

        if ( is_migrate( ptc.q, bounds ) ) {
          //TODO need to encode species into ptc state, otherwise this information is lost
          migrate_buffer.push_back(ptc);
          ptc.set<particle::flag::empty>();
        }

      }



      if constexpr ( sp.is_charged ) {
        WJ *= charge * -grid.delta / dt;

        // check here if needs to output species J
        if ( is_dataexport ) {
          data_out.WJ[sp] = WJ;
          todo::getJ_from_WJ(data_out.J[sp], data_out.WJ[sp] );
        }
      }

      // TODO sort particle array periodically

    }

    migrate( migrate_buffer, neighbors, bounds, comm );
    // TODO load particles back in

    // do this before J is actually used
    auto& J;
    todo::getJ_from_WJ(J, WJ);

  }
}
