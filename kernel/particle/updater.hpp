#ifndef  _PARTICLE_UPDATER_HPP_
#define  _PARTICLE_UPDATER_HPP_

#include "particle/species_predef.hpp"
#include "particle/map.hpp"
#include "particle/array.hpp"

#include "manifold/grid.hpp"

#include "random/rng.hpp"

namespace field {
  template < typename, int, int > struct Field;
}

namespace particle {
  template < int DGrid,
             typename Real,
             template < typename > class PtcSpecs,
             typename ShapeF,
             typename RealJ,
             typename Metric >
  class Updater {
  private:
    const mani::Grid< Real, DGrid >& _localgrid;
    util::Rng<Real> _rng;

    void update_species( species sp,
                         array<Real,PtcSpecs>& sp_ptcs,
                         field::Field<RealJ,3,DGrid>& J,
                         Real dt,
                         const field::Field<Real,3,DGrid>& E,
                         const field::Field<Real,3,DGrid>& B
                         );

  public:
    Updater( const mani::Grid< Real, DGrid >& localgrid, const util::Rng<Real>& rng );

    void operator() ( map<array<Real,PtcSpecs>>& particles,
                      field::Field<RealJ,3,DGrid>& J,
                      const field::Field<Real,3,DGrid>& E,
                      const field::Field<Real,3,DGrid>& B,
                      Real dt, int timestep );
  };

}

#endif
