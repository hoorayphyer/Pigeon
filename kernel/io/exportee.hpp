#ifndef _IO_EXPORTEE_HPP_
#define _IO_EXPORTEE_HPP_

#include "particle/properties.hpp"

namespace mpi { struct CartComm; }

namespace io {
  template < typename RealDS,
             int DGrid,
             typename Real,
             typename ShapeF,
             typename RealJ >
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
    virtual ~FieldExportee() = default;
  };

  template < typename RealDS,
             int DGrid,
             typename Real,
             typename ShapeF,
             typename RealJ >
  struct FieldAction : public FieldExportee<RealDS, DGrid, Real, ShapeF, RealJ> {
    // NOTE we use Real here so RealJ may be downcast. But it is fine since RealJ is mainly to avoid losing precision in current deposition. Here Real is sufficient.
  private:
    using F = apt::array<Real,3> (*) ( apt::Index<DGrid> I,
                                       const mani::Grid<Real,DGrid>& grid,
                                       const field::Field<Real, 3, DGrid>& E,
                                       const field::Field<Real, 3, DGrid>& B,
                                       const field::Field<RealJ, 3, DGrid>& J);
    using H = void (*) ( field::Field<RealDS,3,DGrid>& fds, const mani::Grid<RealDS,DGrid>& grid_ds, int num_comps, const mpi::CartComm& cart );

  public:
    FieldAction( std::string a, int b, F c, H d )
      : varname(a), num_comps(b), impl(c), post_hook(d) {}

    std::string varname {};
    int num_comps = 1;
    F impl = nullptr;
    H post_hook = nullptr;

    virtual std::tuple< std::string, int, field::Field<RealDS,3,DGrid> >
    action( const int ds_ratio,
            const std::optional<mpi::CartComm>& cart_opt,
            const mani::Grid<Real,DGrid>& grid, // local grid
            const mani::Grid<RealDS,DGrid>& grid_ds, // grid of downsampled field
            const int guard_ds,
            const field::Field<Real, 3, DGrid>& E,
            const field::Field<Real, 3, DGrid>& B,
            const field::Field<RealJ, 3, DGrid>& J// J is Jmesh on a replica // TODO double check
            ) override;
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
    virtual ~PtcExportee() = default;
  };

  template < typename RealDS,
             int DGrid,
             typename Real,
             template < typename > class S,
             typename ShapeF
             >
  struct PtcAction : public PtcExportee < RealDS, DGrid, Real, S, ShapeF > {
    using F = apt::array<Real,3> (*) ( const particle::Properties& prop, const typename particle::array<Real,S>::const_particle_type& ptc);
    using H = void (*) ( field::Field<RealDS,3,DGrid>& fds, const mani::Grid<RealDS,DGrid>& grid_ds, int num_comps );

    PtcAction( std::string a, int b, F c, H d )
      : varname(a), num_comps(b), impl(c), post_hook(d) {}

    std::string varname {};
    int num_comps = 1;
    F impl = nullptr;
    H post_hook = nullptr; // here to apply boundary conditions

    virtual std::tuple< std::string, int, field::Field<RealDS,3,DGrid> >
    action( const int ds_ratio,
            const mani::Grid<Real,DGrid>& grid, // local grid
            const mani::Grid<RealDS,DGrid>& grid_ds, // grid of downsampled field
            const int guard_ds,
            particle::species sp,
            const particle::array<Real,S>& ptcs
            ) override;
  };
}


#endif
