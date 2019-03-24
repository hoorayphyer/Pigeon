#include "particle_updater.hpp"
#include "ensemble/ensemble.hpp"

#include "particle/pusher.hpp"
#include "particle/pair_producer.hpp"
#include "particle_properties.hpp"
#include "particle/migration.hpp"

#include "field/field_shape_interplay.hpp"
#include <algorithm>

#include "dynamic_variables.hpp"
#include "parameters.hpp"
#include "kernel/coordinate.hpp"

#include "kernel/grid.hpp"
#include "kernel/shapef.hpp"

namespace particle {
  map<Properties> properties;
}

namespace aperture {
  template < typename Real, int DGrid, int DPtc, typename state_t, typename ShapeF,
             particle::PairScheme pair_scheme, knl::coordsys CS >
  template < bool IsCharged, bool IsRadiative >
  void ParticleUpdater< Real, DGrid, DPtc, state_t, ShapeF, pair_scheme, CS >
  ::update_species( particle::array<Real,3,state_t>& sp_ptcs,
                    particle::map<particle::array<Real,3,state_t>>& particles,
                    Real dt, Real unit_e,
                    const particle::Properties& prop,
                    const field::Field<Real,3,DGrid>& E,
                    const field::Field<Real,3,DGrid>& B,
                    const apt::array< apt::pair<Real>, DGrid >& borders ) {
    if ( sp_ptcs.size() == 0 ) return;

    using namespace particle;
    auto shapef = ShapeF();
    auto charge_over_dt = prop.charge_x * unit_e / dt;

    for ( auto& ptc : sp_ptcs ) { // TODOL check ptc is proxy
      if( ptc.is(flag::empty) ) continue;

      // TODOL should also check IsMassive
      if constexpr ( IsCharged ) {
          auto E_itpl = field::interpolate( E, ptc.q, shapef );
          auto B_itpl = field::interpolate( B, ptc.q, shapef );
          auto&& dp = update_p( ptc, dt, prop.mass_x, E_itpl, B_itpl );

          // TODOL make sure the newly added particles will not participate in the loop
          if constexpr ( IsRadiative ) {
              // TODOL wrap into a factory of pair producer
              if constexpr ( pair_scheme != PairScheme::Disabled ) {
                  auto Rc = calc_Rc( dt, ptc.p, std::move(dp) );
                  auto gamma = std::sqrt( (!prop.mass_x) + apt::sqabs(ptc.p) );
                  if ( is_productive_lepton( ptc, gamma, Rc, _rng ) ) {
                    if constexpr ( pair_scheme == PairScheme::Photon ) {
                        produce_photons( std::back_inserter( particles[species::photon] ),
                                         ptc, gamma, Rc );
                      } else if ( pair_scheme == PairScheme::Instant ) {
                      instant_produce_pairs( std::back_inserter( particles[species::electron] ),
                                             std::back_inserter( particles[species::positron] ),
                                             ptc, gamma, Rc );
                    }
                  }
                }
            }

        }
      else if ( pair_scheme == PairScheme::Photon ) {
        // TODO check is photon
        photon_produce_pairs( std::back_inserter( particles[species::electron] ),
                              std::back_inserter( particles[species::positron] ),
                              ptc );
      }

      // NOTE q is updated, starting from here, particles may be in the guard cells.
      auto&& dq = update_q<CS>( ptc, dt, prop.mass_x );
      // TODOL pusher handle boundary condition. Is it needed?
      if constexpr ( IsCharged )
                     _dJ.deposit( charge_over_dt, ptc.q()-std::move(dq), ptc.q(), shapef ); // TODOL check the 2nd argument

      if ( is_migrate( ptc.q(), borders ) )
        _migrators.emplace_back(std::move(ptc));
    }
  }
}

namespace aperture {

  template < typename Real, int DGrid, int DPtc, typename state_t,
             typename ShapeF,
             particle::PairScheme pair_scheme, knl::coordsys CS
             >
  void ParticleUpdater< Real, DGrid, DPtc, state_t, ShapeF, pair_scheme, CS >
  ::operator() ( field::Field<Real,3,DGrid>& J,
                 particle::map<particle::array<Real,3,state_t>>& particles,
                 const field::Field<Real,3,DGrid>& E,
                 const field::Field<Real,3,DGrid>& B,
                 Real dt, Real unit_e, int timestep ) {
    using namespace particle;
    // borders tell which particles will be migrated. These are simply the boundaries of the bulk grid
    apt::array< apt::pair<Real>, DGrid > borders;
    apt::foreach<0,DGrid>
      ( []( auto& b, const auto& g ) {
          b[LFT] = g.lower();
          b[RGT] = g.upper();
        }, borders, E.mesh().bulk() );

    _dJ.reset();

    // TODO NOTE injection and annihilation are more like free sources and sinks of particles users control, in other words they fit the notion of particle boundary conditions. Do them outside of Updater. In contrast, pair creation is a physical process we model, so it is done inside the updater.
    // inject_particles(); // TODO only coordinate space current is needed in implementing current regulated injection // TODO compatible with ensemble??
    // TODOL annihilation will affect deposition // NOTE one can deposit in the end
    // annihilate_mark_pairs( );

    for ( auto&[ sp, ptcs ] : particles ) {
      const particle::Properties& prop = properties.at(sp);
      auto[ m_x, q_x, is_rad ] = prop;
      // if ( q_x && is_rad )
      //   update_species<true,true>( ptcs, particles, dt, unit_e, E, B, borders );
      // else if ( q_x && !is_rad )
      //   update_species<true,false>( ptcs, particles, dt, unit_e, E, B, borders );
      // else
      //   update_species<false,false>( ptcs, particles, dt, unit_e, E, B, borders );
    }

    migrate( _migrators, _ensemble.inter, borders, timestep );
    for ( auto&& ptc : _migrators ) {
      particles[ptc.template get<species>()].push_back( std::move(ptc) ); // TODOL check this
    }
    _migrators.resize(0);

    _dJ.reduce( _ensemble.chief, _ensemble.intra );

    if ( _cart )
      _dJ.integrate( *_cart );

    // J = getJfromJmesh(_dJ); TODO

  }

}

#include "traits.hpp"
using namespace traits;
namespace aperture {
  template class ParticleUpdater< real_t, DGrid, DPtc, ptc_state_t, knl::shapef_t<shape>,
                                  pair_produce_scheme, coordinate_system>;
}
