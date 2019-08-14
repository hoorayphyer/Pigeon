#ifndef _FIELD_CARTESIAN_UPDATER_HPP_
#define _FIELD_CARTESIAN_UPDATER_HPP_

#include "field/field.hpp"
#include "manifold/grid.hpp"

namespace mpi { struct CartComm; }

namespace field {
  template < typename Real, int DGrid, typename RealJ >
  struct CartesianUpdater {
    void operator() ( Field<Real,3,DGrid>& E,
                      Field<Real,3,DGrid>& B,
                      const Field<RealJ,3,DGrid>& Jmesh,
                      const mani::Grid<Real,DGrid>& grid,
                      const mpi::CartComm& cart,
                      Real dt,
                      Real preJ_factor
                      ) const;
  };
}

#endif
