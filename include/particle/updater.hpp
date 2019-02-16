#ifndef  _PARTICLE_UPDATER_HPP_
#define  _PARTICLE_UPDATER_HPP_

#include "kernel/shape_predef.hpp"
#include "kernel/coordsys_predef.hpp"
#include "particle/pair_produce_predef.hpp"

#include "field/field.hpp"
#include "particle/migration.hpp" // for cParticle
#include "utility/rng.hpp"
#include "parallel/mpi++.hpp"

namespace knl {
  template < int, typename > struct Grid;
}

template < typename, int, int, typename > struct DynamicVars;
template < typename, int > struct Params;

namespace particle {
  // TODO this template parameter list is ugly. Maybe use Policy, for injection, for pair_creation?
  template < typename Real, int DGrid, int DPtc, typename state_t,
             knl::shape Shape, PairScheme pair_scheme,
             knl::coordsys CS, species posion // posion = positron || ion in injection
             >
  class Updater {
  private:
    using WJ_t = Real;
    field::Field< WJ_t, 3, DGrid > _WJ;
    std::vector<cParticle<Real, DPtc, state_t>> _migrators;
    util::Rng<Real> _rng;
    std::array< std::array<std::optional<mpi::InterComm>,2>, DGrid > _intercomms;

    template < species sp, typename Dynvar_t, typename Params_t, typename Grid >
    void update_species( Dynvar_t& dvars, const Params_t& params, const Grid& grid );

    template < int I = 0, typename... Args >
    inline void iterate_species( Args&... args ) {
      if constexpr ( I == 4 ) return; // TODO loop over all species. Would void_t help?
      update_species<static_cast<species>(I)>( args... );
      return iterate_species<I+1>( args... );
    }

  public:
    void operator() ( DynamicVars<Real, DGrid, DPtc, state_t>& dvars,
                      const Params<Real, DGrid>& params,
                      const knl::Grid<DGrid, Real>& grid,
                      const mpi::Comm& ensemble );
  };

}

#endif
