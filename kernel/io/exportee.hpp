#ifndef _IO_EXPORTEE_HPP_
#define _IO_EXPORTEE_HPP_

#include "particle/properties.hpp"

namespace mpi { struct CartComm; }

namespace io {
  template < typename RDS, int DGrid, typename R, typename RJ >
  struct FieldExportee {
    virtual std::tuple< std::string, int, field::Field<RDS,3,DGrid> >
    action ( const int ds_ratio,
             const std::optional<mpi::CartComm>& cart_opt, // in case some global operation is needed
             const mani::Grid<R,DGrid>& grid, // local grid
             const mani::Grid<RDS,DGrid>& grid_ds, // grid of downsampled field
             const int guard_ds,
             const field::Field<R, 3, DGrid>& E,
             const field::Field<R, 3, DGrid>& B,
             const field::Field<RJ, 3, DGrid>& J// J is Jmesh on a primary after reduction
             ) { return {}; }
    virtual ~FieldExportee() = default;
  };

}

namespace io {
  template < typename RDS, int DGrid, typename R, template < typename > class S>
  struct PtcExportee {
    virtual std::tuple< std::string, int, field::Field<RDS,3,DGrid> >
    action ( const int ds_ratio,
             const mani::Grid<R,DGrid>& grid, // local grid
             const mani::Grid<RDS,DGrid>& grid_ds, // grid of downsampled field
             const int guard_ds,
             particle::species sp,
             const particle::array<R,S>& ptcs
             ) { return {}; }
    virtual ~PtcExportee() = default;
  };

}


#endif
