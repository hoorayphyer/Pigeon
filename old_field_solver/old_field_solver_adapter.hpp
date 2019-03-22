#ifndef _OLD_FIELDUPDATER_ADAPTER_HPP_
#define _OLD_FIELDUPDATER_ADAPTER_HPP_

#include "../src/aperture/parameters.hpp"
#include "field/field.hpp"
#include <optional>

namespace mpi { struct CartComm; }

namespace ofs {
  // only applicable in 2D log spherical
  template < int DGrid = 2 >
  struct OldFieldUpdater {
    static_assert( DGrid == 2 );
  private:
    const mpi::CartComm& _cart;

  public:
    using field_type = field::Field<double,3,DGrid>;

    OldFieldUpdater( const Params<double>& params,
                     const mpi::CartComm& cart,
                     const knl::Grid<double,DGrid>& local_grid,
                     apt::array< apt::pair<bool>, DGrid > is_at_boundary,
                     apt::array< apt::pair<std::optional<int>>, DGrid > neigh_cart_ranks,
                     int guard
                     );

    void operator() ( field_type& E,
                      field_type& B,
                      const field_type& J,
                      typename field_type::element_type dt, int timestep );
  };
}

#endif
