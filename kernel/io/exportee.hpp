#ifndef _IO_EXPORTEE_HPP_
#define _IO_EXPORTEE_HPP_

namespace io {
  template < typename RealDS,
             int DGrid,
             typename Real,
             typename ShapeF,
             typename RealJ,
             typename Metric >
  struct FieldExportee {
    virtual std::tuple< std::string, int, field::Field<RealDS,3,DGrid> >
    action ( const int ds_ratio,
             const std::optional<mpi::CartComm>& cart_opt, // in case some global operation is needed
             const mani::Grid<Real,DGrid>& grid, // local grid
             const mani::Grid<RealDS,DGrid>& grid_ds, // grid of downsampled field
             const int guard_ds,
             const field::Field<Real, 3, DGrid>& E,
             const field::Field<Real, 3, DGrid>& B,
             const field::Field<RealJ, 3, DGrid>& J// J is Jmesh on a primary after reduction
             ) { return {}; }
  };
}

namespace io {
  template < typename RealDS,
             int DGrid,
             typename Real,
             template < typename > class S,
             typename ShapeF
             >
  struct PtcExportee {
    virtual std::tuple< std::string, int, field::Field<RealDS,3,DGrid> >
    action ( const int ds_ratio,
             const mani::Grid<Real,DGrid>& grid, // local grid
             const mani::Grid<RealDS,DGrid>& grid_ds, // grid of downsampled field
             const int guard_ds,
             particle::species sp,
             const particle::array<Real,S>& ptcs
             ) { return {}; }
  };
}


#endif
