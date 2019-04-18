#ifndef  _PARTICLE_UPDATER_HPP_
#define  _PARTICLE_UPDATER_HPP_

#include "kernel/coordsys_predef.hpp"
#include "particle/species_predef.hpp"
#include "particle/map.hpp"
#include "particle/array.hpp"

#include "kernel/grid.hpp"

#include "utility/rng.hpp"

namespace field {
  template < typename, int, int > struct Field;
}

namespace particle {
  template < typename Real, int DGrid, int DPtc, typename state_t, typename ShapeF,
             typename RealJ, knl::coordsys CS >
  class ParticleUpdater {
  private:
    const knl::Grid< Real, DGrid >& _localgrid;
    util::Rng<Real> _rng;

    void update_species( species sp,
                         array<Real,3,state_t>& sp_ptcs,
                         field::Field<RealJ,3,DGrid>& J,
                         Real dt,
                         const field::Field<Real,3,DGrid>& E,
                         const field::Field<Real,3,DGrid>& B
                         );

  public:
    ParticleUpdater( const knl::Grid< Real, DGrid >& localgrid, const util::Rng<Real>& rng );

    void operator() ( map<array<Real,3,state_t>>& particles,
                      field::Field<RealJ,3,DGrid>& J,
                      const field::Field<Real,3,DGrid>& E,
                      const field::Field<Real,3,DGrid>& B,
                      Real dt, int timestep );
  };

}

#endif
