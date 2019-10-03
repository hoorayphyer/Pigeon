#ifndef _FIELD_ACTION_HPP_
#define _FIELD_ACTION_HPP_

#include "apt/action_base.hpp"
#include "field/field.hpp"
#include "manifold/grid.hpp"

namespace mpi { struct CartComm; }

namespace field {
  template < typename Real, int DGrid, typename RealJ >
  struct Action : public apt::ActionBase<DGrid> {
  private:
    bool _orig_EB = true; // this will tell the simulator to restore values in overlapping areas before applying this action

  public:
    Action& require_original_EB(bool x) { _orig_EB = x; return *this; }

    const auto& orig_EB() const noexcept { return _orig_EB; }

    virtual ~Action() {};

    virtual Action* Clone() const = 0; // covariant return types, see Modern C++ Design

    virtual void operator() ( Field<Real,3,DGrid>& E,
                              Field<Real,3,DGrid>& B,
                              const Field<RealJ,3,DGrid>& Jmesh,
                              const mani::Grid<Real,DGrid>& grid,
                              const mpi::CartComm& cart,
                              int timestep,
                              Real dt
                              ) const = 0;
  };
}

#endif
