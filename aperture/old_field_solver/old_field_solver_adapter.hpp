#ifndef _OLD_FIELDUPDATER_ADAPTER_HPP_
#define _OLD_FIELDUPDATER_ADAPTER_HPP_

#include "field/field.hpp"
#include "kernel/grid.hpp"
#include "../abstract_field_updater.hpp"
#include <optional>

namespace mpi { struct CartComm; }

namespace ofs {
  // only applicable in 2D log spherical
  template < int DGrid = 2 >
  struct OldFieldUpdater : aperture::AbstractFieldUpdater<double,DGrid> {
    static_assert( DGrid == 2 );
  private:
    const mpi::CartComm& _cart;

  public:
    using field_type = field::Field<double,3,DGrid>;


    // apt::array< apt::pair<std::optional<int>>, DGrid > neigh_cart_ranks; // TODOL currently used in old_field_solver, and link_neighbors
    OldFieldUpdater( const mpi::CartComm& cart,
                     const knl::Grid<double,DGrid>& local_grid,
                     apt::array< apt::pair<bool>, DGrid > is_at_boundary,
                     int guard );

    virtual void operator() ( field_type& E,
                              field_type& B,
                              const field_type& J,
                              double dt, int timestep ) override;
  };
}

#endif
