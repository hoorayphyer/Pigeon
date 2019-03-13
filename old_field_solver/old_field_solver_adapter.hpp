#ifndef _OLD_FIELDUPDATER_ADAPTER_HPP_
#define _OLD_FIELDUPDATER_ADAPTER_HPP_

#include "parameters.hpp"
#include "field/field.hpp"
#include <mpi.h>

namespace parallel { template < int > struct Locale; }
namespace mpi { struct Comm; }

namespace ofs {
  // only applicable in 2D log spherical
  struct OldFieldUpdater {
  private:
    const mpi::Comm& _comm;

  public:
    static constexpr int DGrid = 2;
    using field_type = field::Field<double,3,DGrid>;

    OldFieldUpdater( const Params<double>& params,
                     const mpi::Comm& comm,
                     const knl::Grid<double,DGrid>& local_grid,
                     const parallel::Locale<DGrid>& locale,
                     int guard
                     );

    // TODO factor of 4\pi on J ?
    void operator() ( field_type& E,
                      field_type& B,
                      const field_type& J,
                      field_type::element_type dt, int timestep );
  };
}

#endif
