#ifndef  _APERTURE_ABSTRACT_PARTICLE_UPDATER_HPP_
#define  _APERTURE_ABSTRACT_PARTICLE_UPDATER_HPP_

#include "field/field.hpp"
#include "particle/array.hpp"
#include "particle/map.hpp"

namespace aperture {
  template < typename Real, int DGrid, typename state_t >
  struct AbstractParticleUpdater {
    virtual void operator() ( field::Field<Real,3,DGrid>& J,
                              particle::map<particle::array<Real,3,state_t>>& particles,
                              const field::Field<Real,3,DGrid>& E,
                              const field::Field<Real,3,DGrid>& B,
                              Real dt,Real unit_e, int timestep ) = 0;
  };
}

#endif
