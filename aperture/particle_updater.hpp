#ifndef  _PARTICLE_UPDATER_HPP_
#define  _PARTICLE_UPDATER_HPP_

#include "abstract_particle_updater.hpp"

#include "kernel/coordsys_predef.hpp"
#include "particle/species_predef.hpp"
#include "particle/pair_produce_predef.hpp"
#include "particle_properties.hpp"
#include "particle/map.hpp"
#include "particle/array.hpp"

#include "kernel/grid.hpp"
#include "particle/c_particle.hpp"
#include "field/mesh_shape_interplay.hpp"

#include "utility/rng.hpp"
#include "parallel/mpi++.hpp"

namespace dye {
  template < int > struct Ensemble;
}

namespace particle {
  // TODO this template parameter list is ugly. Maybe use Policy, for injection, for pair_creation?
  template < typename Real, int DGrid, int DPtc, typename state_t, typename ShapeF,
             typename Real_dJ,
             PairScheme pair_scheme, knl::coordsys CS >
  class ParticleUpdater : public aperture::AbstractParticleUpdater<Real, DGrid, state_t>{
  private:
    const knl::Grid< Real, DGrid >& _localgrid;
    field::Standard_dJ_Field< Real_dJ, 3, DGrid, ShapeF > _dJ;
    std::vector<cParticle<Real, DPtc, state_t>> _migrators;
    util::Rng<Real> _rng;
    const std::optional<mpi::CartComm>& _cart;
    const dye::Ensemble<DGrid>& _ensemble;

    template < bool IsCharged >
    void update_species( array<Real,3,state_t>& sp_ptcs,
                         Real dt, Real unit_e,
                         const Properties& prop,
                         const field::Field<Real,3,DGrid>& E,
                         const field::Field<Real,3,DGrid>& B,
                         const apt::array< apt::pair<Real>, DGrid >& borders
                         );

  public:
    ParticleUpdater( const knl::Grid< Real, DGrid >& localgrid, const util::Rng<Real>& rng, const std::optional<mpi::CartComm>& cart, const dye::Ensemble<DGrid>& ensemble );

    virtual void operator() ( field::Field<Real,3,DGrid>& J,
                              map<array<Real,3,state_t>>& particles,
                              const field::Field<Real,3,DGrid>& E,
                              const field::Field<Real,3,DGrid>& B,
                              const apt::array< apt::pair<Real>, DGrid >& borders,
                              Real dt,Real unit_e, int timestep ) override;
  };

}

#include "traits.hpp"
#include "kernel/shapef.hpp"
namespace aperture {
  using namespace traits;
  template < typename Real, int DGrid, typename state_t >
  using ParticleUpdater = particle::ParticleUpdater<Real, DGrid, 3, state_t, knl::shapef_t<shape>, real_dj_t, pair_produce_scheme, coordinate_system >;
}

#endif
