#ifndef  _PARTICLE_UPDATER_HPP_
#define  _PARTICLE_UPDATER_HPP_

#include "kernel/shapef.hpp"
#include "kernel/coordsys_predef.hpp"
#include "particle/species_predef.hpp"
#include "particle/pair_produce_predef.hpp"

#include "particle/c_particle.hpp"
#include "field/field_shape_interplay.hpp"

#include "utility/rng.hpp"
#include "parallel/mpi++.hpp"

#include "dynamic_variables.hpp"
#include "parameters.hpp"

namespace aperture {
  template < int > struct Ensemble;

  // TODO this template parameter list is ugly. Maybe use Policy, for injection, for pair_creation?
  template < typename Real, int DGrid, int DPtc, typename state_t,
             knl::shape Shape, particle::PairScheme pair_scheme,
             knl::coordsys CS
             >
  class ParticleUpdater {
  public:
    using DynaVars_t = DynamicVars<Real, DGrid, DPtc, state_t>;
    using Params_t = Params<Real>;

  private:
    using ShapeF = knl::shapef_t<Shape>;

    field::dJ_Field< long double, 3, DGrid > _dJ; // TODO long double is hard coded
    std::vector<particle::cParticle<Real, DPtc, state_t>> _migrators;
    util::Rng<Real> _rng;
    apt::array< apt::pair<std::optional<mpi::InterComm>>, DGrid > _intercomms;

  public:
    // TODO constructor?


    void operator() ( int timestep, DynaVars_t& dvars, const Params_t& params,
                      const std::optional<mpi::CartComm>& cart,
                      const Ensemble<DGrid>& ensemble );
  };

}

#endif
