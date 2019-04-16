#ifndef  _PARTICLE_UPDATER_HPP_
#define  _PARTICLE_UPDATER_HPP_

#include "kernel/coordsys_predef.hpp"
#include "particle/species_predef.hpp"
#include "particle/map.hpp"
#include "particle/array.hpp"

#include "kernel/grid.hpp"
#include "field/current_deposition.hpp"

#include "utility/rng.hpp"

namespace particle {
  template < typename Real, int DGrid, int DPtc, typename state_t, typename ShapeF,
             typename Real_dJ, knl::coordsys CS >
  class ParticleUpdater {
  private:
    const knl::Grid< Real, DGrid >& _localgrid;
    field::Standard_dJ_Field< Real_dJ, 3, DGrid, ShapeF > _dJ;
    util::Rng<Real> _rng;

    using ReturnType = decltype(_dJ.integrate());

    void update_species( species sp,
                         array<Real,3,state_t>& sp_ptcs,
                         Real dt, Real unit_e,
                         const field::Field<Real,3,DGrid>& E,
                         const field::Field<Real,3,DGrid>& B
                         );

  public:
    ParticleUpdater( const knl::Grid< Real, DGrid >& localgrid, const util::Rng<Real>& rng );

    ReturnType operator() ( map<array<Real,3,state_t>>& particles,
                 const field::Field<Real,3,DGrid>& E,
                 const field::Field<Real,3,DGrid>& B,
                 Real dt,Real unit_e, int timestep );
  };

}

#endif
