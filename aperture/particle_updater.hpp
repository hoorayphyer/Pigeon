#ifndef  _PARTICLE_UPDATER_HPP_
#define  _PARTICLE_UPDATER_HPP_

#include "abstract_particle_updater.hpp"

#include "kernel/coordsys_predef.hpp"
#include "particle/species_predef.hpp"
#include "particle/pair_produce_predef.hpp"

#include "particle/c_particle.hpp"
#include "field/field_shape_interplay.hpp"

#include "utility/rng.hpp"
#include "parallel/mpi++.hpp"

#include "dynamic_variables.hpp"
#include "parameters.hpp"

namespace particle { struct Properties; }

namespace aperture {
  template < int > struct Ensemble;

  // TODO this template parameter list is ugly. Maybe use Policy, for injection, for pair_creation?
  template < typename Real, int DGrid, int DPtc, typename state_t, typename ShapeF,
             particle::PairScheme pair_scheme, knl::coordsys CS >
  class ParticleUpdater : public AbstractParticleUpdater<Real, DGrid, state_t>{
  private:
    field::dJ_Field< long double, 3, DGrid > _dJ; // TODO long double is hard coded
    std::vector<particle::cParticle<Real, DPtc, state_t>> _migrators;
    util::Rng<Real> _rng;
    const std::optional<mpi::CartComm>& _cart;
    const Ensemble<DGrid>& _ensemble;

    template < bool IsCharged, bool IsRadiative >
    void update_species( particle::array<Real,3,state_t>& sp_ptcs,
                         particle::map<particle::array<Real,3,state_t>>& particles,
                         Real dt, Real unit_e,
                         const particle::Properties& prop,
                         const field::Field<Real,3,DGrid>& E,
                         const field::Field<Real,3,DGrid>& B,
                         const apt::array< apt::pair<Real>, DGrid >& borders );




  public:
    ParticleUpdater( const knl::Grid< Real, DGrid, knl::grid1d::Clip >& localgrid, const util::Rng<Real>& rng, const std::optional<mpi::CartComm>& cart, const Ensemble<DGrid>& ensemble );

    virtual void operator() ( field::Field<Real,3,DGrid>& J,
                              particle::map<particle::array<Real,3,state_t>>& particles,
                              const field::Field<Real,3,DGrid>& E,
                              const field::Field<Real,3,DGrid>& B,
                              Real dt,Real unit_e, int timestep ) override;
  };

}

#endif
