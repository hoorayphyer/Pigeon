#ifndef  _PARTICLE_UPDATER_HPP_
#define  _PARTICLE_UPDATER_HPP_

#include "kernel/shapef.hpp"
#include "kernel/coordsys_predef.hpp"
#include "particle/species_predef.hpp"
#include "particle/pair_produce_predef.hpp"

#include "particle/migration.hpp"
#include "utility/rng.hpp"
#include "parallel/mpi++.hpp"

#include "dynamic_variables.hpp"
#include "parameters.hpp"

namespace particle {
  // TODO this template parameter list is ugly. Maybe use Policy, for injection, for pair_creation?
  template < typename Real, int DGrid, int DPtc, typename state_t,
             knl::shape Shape, PairScheme pair_scheme,
             knl::coordsys CS, species posion // posion = positron || ion in injection
             >
  class Updater {
  public:
    using DynaVars_t = DynamicVars<Real, DGrid, DPtc, state_t>;
    using Params_t = Params<Real, DGrid>;

  private:
    using WJ_t = long double; // TODO
    using ShapeF = knl::shapef_t<Shape>;

    field::Field< WJ_t, 3, DGrid > _WJ;
    std::vector<cParticle<Real, DPtc, state_t>> _migrators;
    util::Rng<Real> _rng;
    apt::array< apt::pair<std::optional<mpi::InterComm>>, DGrid > _intercomms;

  public:
    void operator() ( int timestep, DynaVars_t& dvars, const Params_t& params,
                      const std::optional<mpi::Comm>& ensemble );
  };

}

#endif
