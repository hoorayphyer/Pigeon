#ifndef _FIELDUPDATER_HPP_
#define _FIELDUPDATER_HPP_

#include "field/field.hpp"
#include "manifold/grid.hpp"
#include <optional>

namespace mpi { struct CartComm; }

namespace field {
  template < typename Real, int DGrid, typename RealJ >
  struct Updater {
  private:
    const mpi::CartComm& _cart;

  public:
    using field_type = field::Field<double,3,DGrid>;
    using J_type = field::Field<RealJ,3,DGrid>;

    // apt::array< apt::pair<std::optional<int>>, DGrid > neigh_cart_ranks; // TODOL currently used in old_field_solver, and link_neighbors
    Updater( const mpi::CartComm& cart,
             const mani::Grid<double,DGrid>& local_grid,
             apt::array< apt::pair<bool>, DGrid > is_at_boundary,
             int guard );

    void operator() ( field_type& E,
                      field_type& B,
                      const J_type& Jmesh,
                      double dt, int timestep );
  };
}

#endif
