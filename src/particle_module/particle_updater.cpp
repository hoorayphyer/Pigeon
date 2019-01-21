#include "particle_module/particle_updater.hpp"
#inlcude "particle_module/pusher.hpp"
#include "particle_module/depositer.hpp"

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



  template < sf::shape S >
  void update( ) {
    auto& WJ;
    WJ.zero_out();
    EnsBroadcastFields( Efield, Bfield );

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
          auto&& dp = update_p();

          // TODO make sure the newly added particles will not participate in the loop
          if constexpr ( sp.is_radiative ) {
            if constexpr ( photon == pair ) {
              procude_photon( todo::find_Rc( std::move(dp) ), photons );
            } else if ( instant == pair ) {
              electrons, positrons;
            }
          }

          auto&& dq = update_q(); // NOTE q is updated
          // pusher handle boundary condition

          depositWJ<S>( WJ, ptc, std::move(dq), grid );

        } else if ( sp == photon ) {
          update_q(); // NOTE q is updated
          photon.producepairs;

        } else {
          static_assert( sp_not_here );
        }

        sort_into_communication_buffer(ptc);


      }
      //-------------------------------
      communicate_particles_for_this_species( particles );



      if constexpr ( sp.is_charged ) {
        WJ *= charge * -grid.delta / dt;

        // check here if needs to output species J
        if ( is_dataexport ) {
          data_out.WJ[sp] = WJ;
          todo::getJ_from_WJ(data_out.J[sp], data_out.WJ[sp] );
        }
      }

    }

    // do this before J is actually used
    auto& J;
    todo::getJ_from_WJ(J, WJ);

  }
}
