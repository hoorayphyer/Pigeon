#ifndef _IO_EXPORTEE_BY_FUNCTION_HPP_
#define _IO_EXPORTEE_BY_FUNCTION_HPP_

#include "io/exportee.hpp"
#include "msh/mesh_shape_interplay.hpp"
#include <cassert>

namespace io {
  constexpr int POW( int B, int E ) {
    if ( E == 0 ) return 1;
    else return B * POW(B,E-1);
  };
}

namespace io {
  template < typename RealDS,
             int DGrid,
             typename Real,
             typename RealJ >
  struct FexpTbyFunction : public FieldExportee<RealDS, DGrid, Real, RealJ> {
    // NOTE we use Real here so RealJ may be downcast. But it is fine since RealJ is mainly to avoid losing precision in current deposition. Here Real is sufficient.
    using F = apt::array<Real,3> (*) ( apt::Index<DGrid> I,
                                       const mani::Grid<Real,DGrid>& grid,
                                       const field::Field<Real, 3, DGrid>& E,
                                       const field::Field<Real, 3, DGrid>& B,
                                       const field::Field<RealJ, 3, DGrid>& J);
    using H = void (*) ( field::Field<RealDS,3,DGrid>& fds, const mani::Grid<RealDS,DGrid>& grid_ds, int num_comps, const mpi::CartComm& cart );

    FexpTbyFunction( std::string a, int b, F c, H d )
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
            ) override {
      if ( num_comps < 1 ) return {};
      assert( impl != nullptr );

      field::Field<RealDS,3,DGrid> fds( { mani::dims(grid_ds), guard_ds } );
      for ( int i = 0; i < num_comps; ++i ) {
        for ( int dim = 0; dim < DGrid; ++dim )
          fds.set_offset( i, dim, MIDWAY );
      }

      apt::Index<DGrid> subext;
      for ( int i = 0; i < DGrid; ++i ) subext[i] = ds_ratio;

      for ( const auto& Ids : apt::Block( fds.mesh().bulk_dims() ) ) {
        for ( const auto& Isub : apt::Block(subext) ) {
          apt::Index<DGrid> I;
          for ( int i = 0; i < DGrid; ++i )
            I[i] = ds_ratio * Ids[i] + Isub[i];

          auto val = impl( I, grid, E, B, J );

          for ( int comp = 0; comp < num_comps; ++comp )
            fds[comp](Ids) += val[comp];
        }
      }
      {
        int factor = POW(ds_ratio, DGrid);
        for ( int i = 0; i < num_comps; ++i ) {
          for ( auto& x : fds[i].data() ) x /= factor;
        }
      }

      if ( post_hook != nullptr && cart_opt ) post_hook( fds, grid_ds, num_comps, *cart_opt );

      return std::make_tuple( varname, num_comps, fds );
    }
  };

}

namespace io {
  template < typename RealDS,
             int DGrid,
             typename Real,
             template < typename > class S,
             typename ShapeF
             >
  struct PexpTbyFunction : public PtcExportee < RealDS, DGrid, Real, S > {
    using F = apt::array<Real,3> (*) ( const particle::Properties& prop, const typename particle::array<Real,S>::const_particle_type& ptc);
    using H = void (*) ( field::Field<RealDS,3,DGrid>& fds, const mani::Grid<RealDS,DGrid>& grid_ds, int num_comps );

    PexpTbyFunction( std::string a, int b, F c, H d )
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
            ) override {
      if ( num_comps < 1 ) return {};
      assert( impl != nullptr );

      ShapeF sf;

      field::Field<RealDS,3,DGrid> fds( { mani::dims(grid_ds), guard_ds } );
      for ( int i = 0; i < num_comps; ++i ) {
        for ( int dim = 0; dim < DGrid; ++dim )
          fds.set_offset( i, dim, MIDWAY );
      }

      const auto& prop = particle::properties[sp];
      for ( const auto& ptc : ptcs ) {
        if ( !ptc.is(particle::flag::exist) ) continue;
        msh::deposit( fds, ptc.frac(), impl( prop, ptc ), msh::to_standard(grid_ds, ptc.q()), sf);
      }

      if ( post_hook != nullptr ) post_hook( fds, grid_ds, num_comps );

      return std::make_tuple(varname, num_comps, fds);
    }


  };
}

#endif