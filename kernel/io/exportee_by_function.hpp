#ifndef _IO_EXPORTEE_BY_FUNCTION_HPP_
#define _IO_EXPORTEE_BY_FUNCTION_HPP_

#include "io/exportee.hpp"
// #include "msh/mesh_shape_interplay_impl.hpp" // FIXME need deposit with induced shape
#include <cassert>

// FIXME check if downsampled field can use a better mesh. Range may not work because of the scaling

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
                                       const apt::Grid<Real,DGrid>& grid,
                                       const field::Field<Real, 3, DGrid>& E,
                                       const field::Field<Real, 3, DGrid>& B,
                                       const field::Field<RealJ, 3, DGrid>& J);
    using H = void (*) ( field::Field<RealDS,3,DGrid>& fds, const apt::Grid<RealDS,DGrid>& grid_ds, int num_comps, const mpi::CartComm& cart );

    FexpTbyFunction( std::string a, int b, F c, H d )
      : impl(c), post_hook(d) { this->_name = a; this->_num_comps = b; }

    F impl = nullptr;
    H post_hook = nullptr;

    virtual field::Field<RealDS,3,DGrid>
    action( const int ds_ratio,
            const std::optional<mpi::CartComm>& cart_opt,
            const apt::Grid<Real,DGrid>& grid, // local grid
            const apt::Grid<RealDS,DGrid>& grid_ds, // grid of downsampled field
            const int guard_ds,
            const field::Field<Real, 3, DGrid>& E,
            const field::Field<Real, 3, DGrid>& B,
            const field::Field<RealJ, 3, DGrid>& J// J is Jmesh on a replica // TODO double check
            ) override {
      const auto& num_comps = this->_num_comps;
      if ( num_comps < 1 ) return {};
      assert( impl != nullptr );

      field::Field<RealDS,3,DGrid> fds{ apt::make_range( {}, apt::dims(grid_ds), guard_ds ) };
      for ( int i = 0; i < num_comps; ++i ) {
        for ( int dim = 0; dim < DGrid; ++dim )
          fds.set_offset( i, dim, MIDWAY );
      }

      apt::Index<DGrid> subext;
      for ( int i = 0; i < DGrid; ++i ) subext[i] = ds_ratio;

      for ( const auto& Ids : apt::Block( {}, apt::range::size(fds.mesh().range()) ) ) {
        for ( const auto& Isub : apt::Block({}, subext) ) {
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

      return fds;
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
    using H = void (*) ( field::Field<RealDS,3,DGrid>& fds, const apt::Grid<RealDS,DGrid>& grid_ds, int num_comps );

    PexpTbyFunction( std::string a, int b, F c, H d )
      : impl(c), post_hook(d) {this->_name = a; this->_num_comps = b; }

    F impl = nullptr;
    H post_hook = nullptr; // here to apply boundary conditions

    virtual field::Field<RealDS,3,DGrid>
    action( const int ds_ratio,
            const apt::Grid<Real,DGrid>& grid, // local grid
            const apt::Grid<RealDS,DGrid>& grid_ds, // grid of downsampled field
            const int guard_ds,
            const particle::Properties& prop,
            const particle::array<Real,S>& ptcs
            ) override {
      const auto& num_comps = this->_num_comps;
      if ( num_comps < 1 ) return {};
      assert( impl != nullptr );

      ShapeF sf;

      field::Field<RealDS,3,DGrid> fds{ apt::make_range( {}, apt::dims(grid_ds), guard_ds ) };
      for ( int i = 0; i < num_comps; ++i ) {
        for ( int dim = 0; dim < DGrid; ++dim )
          fds.set_offset( i, dim, MIDWAY );
      }

      for ( const auto& ptc : ptcs ) {
        if ( !ptc.is(particle::flag::exist) ) continue;
        // msh::deposit( fds, ptc.frac(), impl( prop, ptc ), msh::to_standard(grid_ds, ptc.q()), sf);
        // simple nearest-cell policy is used. NOTE to compare with MIDWAY
        auto val = impl(prop,ptc);
        apt::Index<DGrid> I;
        for ( int i = 0; i < DGrid; ++i )
          I[i] = ( ptc.q(i) - grid_ds[i].lower() ) / grid_ds[i].delta();
        for ( int i = 0; i < num_comps; ++i )
          fds[i](I) += val[i] * ptc.frac();

      }

      if ( post_hook != nullptr ) post_hook( fds, grid_ds, num_comps );

      return fds;
    }


  };
}

#endif
