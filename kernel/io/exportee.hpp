#ifndef _IO_EXPORTEE_HPP_
#define _IO_EXPORTEE_HPP_

#include "apt/grid.hpp"
#include "field/field.hpp"
#include "particle/array.hpp"
#include "particle/map.hpp"
#include "particle/properties.hpp"
#include <optional>

namespace mpi { struct CartComm; }

namespace io {
  struct ExporteeBase {
  protected:
    std::string _name;
    int _num_comps = 1;

  public:
    inline const std::string& name() const noexcept { return _name; }
    inline const int& num_comps() const noexcept { return _num_comps; }
  };


  template < typename RDS, int DGrid, typename R, typename RJ >
  struct FieldExportee : public ExporteeBase {
    virtual field::Field<RDS,3,DGrid>
    action ( const int ds_ratio,
             const std::optional<mpi::CartComm>& cart_opt, // in case some global operation is needed
             const apt::Grid<R,DGrid>& grid, // local grid
             const apt::Grid<RDS,DGrid>& grid_ds, // grid of downsampled field
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
  struct PtcExportee : public ExporteeBase {
    virtual field::Field<RDS,3,DGrid>
    action ( const int ds_ratio,
             const apt::Grid<R,DGrid>& grid, // local grid
             const apt::Grid<RDS,DGrid>& grid_ds, // grid of downsampled field
             const int guard_ds,
             const particle::Properties& prop,
             const particle::array<R,S>& ptcs
             ) { return {}; }
    virtual ~PtcExportee() = default;
  };

}


#endif
