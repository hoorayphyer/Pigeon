#ifndef  _PARTICLE_UPDATER_HPP_
#define  _PARTICLE_UPDATER_HPP_

#include "kernel/coordsys_predef.hpp"
#include "particle/species_predef.hpp"
#include "particle/map.hpp"
#include "particle/array.hpp"

#include "kernel/grid.hpp"
#include "particle/c_particle.hpp"
#include "field/current_deposition.hpp"

#include "utility/rng.hpp"
#include "parallel/mpi++.hpp"

namespace dye {
  template < int > struct Ensemble;
}

namespace particle {
  template < typename Real, int DGrid, int DPtc, typename state_t, typename ShapeF,
             typename Real_dJ, knl::coordsys CS >
  class ParticleUpdater {
  private:
    const knl::Grid< Real, DGrid >& _localgrid;
    field::Standard_dJ_Field< Real_dJ, 3, DGrid, ShapeF > _dJ;
    std::vector<cParticle<Real, DPtc, state_t>> _migrators;
    util::Rng<Real> _rng;
    const std::optional<mpi::CartComm>& _cart;
    const dye::Ensemble<DGrid>& _ensemble;

    using ReturnType = decltype(_dJ.integrate());

    void update_species( species sp,
                         array<Real,3,state_t>& sp_ptcs,
                         Real dt, Real unit_e,
                         const field::Field<Real,3,DGrid>& E,
                         const field::Field<Real,3,DGrid>& B,
                         const apt::array< apt::pair<Real>, DGrid >& borders
                         );

  public:
    ParticleUpdater( const knl::Grid< Real, DGrid >& localgrid, const util::Rng<Real>& rng, const std::optional<mpi::CartComm>& cart, const dye::Ensemble<DGrid>& ensemble );

    ReturnType operator() ( map<array<Real,3,state_t>>& particles,
                 const field::Field<Real,3,DGrid>& E,
                 const field::Field<Real,3,DGrid>& B,
                 const apt::array< apt::pair<Real>, DGrid >& borders,
                 Real dt,Real unit_e, int timestep );
  };

}

#endif
