#include "field/updater.hpp"

#include "field/params.hpp"
#include "field/field.hpp"
#include "field/sync.hpp"

#include "old_field_solver/FieldUpdater.h"
#include "old_field_solver/FieldCommunicator.h"
#include "old_field_solver/FiniteDiff.h"
#include "old_field_solver/CoordSystem.h"

#include "mpipp/mpi++.hpp"
#include <cmath>
#include <memory>

// NOTE these are modulo max_omega
void Rotating_Monopole_LogSph( FBC& bdry ) {
  bdry.B1 = [] (Scalar r_log, Scalar , Scalar )-> Scalar { return 1.0 / std::exp(2*r_log); };
  bdry.B2 = [] (Scalar, Scalar, Scalar)-> Scalar { return 0.0; };
  bdry.B3 = [] (Scalar, Scalar, Scalar)-> Scalar { return 0.0; };

  // find out E by E = -( Omega x r ) x B
  bdry.E1 = []( Scalar, Scalar, Scalar)-> Scalar { return 0.0; };
  bdry.E2 = []( Scalar r_log, Scalar theta, Scalar )-> Scalar { return - std::exp(-r_log) * std::sin( theta ); };
  bdry.E3 = []( Scalar, Scalar, Scalar )-> Scalar { return 0.0; };
}

void Rotating_Dipole_LogSph( FBC& bdry ) {
  bdry.B1 = [] (Scalar x, Scalar theta, Scalar)-> Scalar { return 2.0 * std::cos(theta) * std::exp(-3*x); };
  bdry.B2 = [] (Scalar x, Scalar theta, Scalar)-> Scalar { return std::sin(theta) * std::exp(-3*x); };
  bdry.B3 = [] (Scalar, Scalar, Scalar)-> Scalar { return 0.0; };

  // find out E by E = -( Omega x r ) x B
  bdry.E1 = []( Scalar r_log, Scalar theta, Scalar )-> Scalar {
              return std::exp(- 2 * r_log) * std::sin( theta ) * std::sin(theta);};
  bdry.E2 = []( Scalar r_log, Scalar theta, Scalar )-> Scalar {
              return - std::exp(- 2 * r_log) * std::sin( 2*theta );};
  bdry.E3 = []( Scalar, Scalar, Scalar )-> Scalar {return 0.0;};
}

void RestoreJToRealSpace( VectorField<Scalar>& JField, const Grid& grid ) {
  // static so that it can be used in lambda functions
  static typename CoordToScales<CoordType::LOG_SPHERICAL>::type coord;

  // define a function pointer.
  Scalar (*h_func) ( std::array<Scalar,3> q ) = nullptr;

  for ( int comp = 0; comp < 3; ++comp ) {
    Index stagJ;
    // NOTE only non-capture lambda can be assigned to function pointer
    stagJ = GetStagProperty(FieldType::ETYPE, comp);
    switch(comp) {
    case 0 : h_func =
        [] ( std::array<Scalar,3> q )-> Scalar {
          return coord.h2( q[0], q[1], q[2] ) * coord.h3( q[0], q[1], q[2] ); }; break;
    case 1 : h_func =
        [] ( std::array<Scalar,3> q )-> Scalar {
          return coord.h3( q[0], q[1], q[2] ) * coord.h1( q[0], q[1], q[2] ); }; break;
    case 2 : h_func =
        [] ( std::array<Scalar,3> q )-> Scalar {
          return coord.h1( q[0], q[1], q[2] ) * coord.h2( q[0], q[1], q[2] ); }; break;
    }

    std::array<Scalar,3> q = { 0, 0, 0 };
    for ( int k = grid.guard[2]; k < grid.dims[2] - grid.guard[2]; ++k ) {
      q[2] = grid.pos( 2, k, stagJ[2] );
      for ( int j = grid.guard[1]; j < grid.dims[1] - grid.guard[1]; ++j ) {
        q[1] = grid.pos( 1, j, stagJ[1] );
        for( int i = grid.guard[0]; i < grid.dims[0] - grid.guard[0]; ++i ) {
          q[0] = grid.pos( 0, i, stagJ[0]);

          if ( std::abs(h_func(q)) > 1e-12 )
            JField( comp, i, j, k ) /= h_func(q);
          else
            JField( comp, i, j, k ) = 0.0;
        }
      }
    }
  }

  return;
}

namespace field {
  VectorField<Scalar> Efield;
  VectorField<Scalar> Bfield;
  VectorField<Scalar> current;

  FUParams fuparams;

  std::unique_ptr<FieldCommunicator> fc{};
  std::unique_ptr<FiniteDiff> fd{};
  std::unique_ptr<FieldUpdater> fu{};

  // Set field BCs for rotating conductor in 2D log spherical coordinates with a damping layer
  template < typename FBCs >
  void Set_FBC_RC_2DLogSph_Damp( FBCs& fieldBC, int guard, Scalar (*omega_t)(Scalar) ) {
    // field BC, specialized for pulsar
    { auto& fbc = fieldBC[LOWER_1];
      fbc.type = FieldBCType::ROTATING_CONDUCTOR;
      fbc.indent = ofs::indent[LOWER_1];
      fbc.ft = omega_t;
      if ( ofs::magnetic_pole == 1 )
        Rotating_Monopole_LogSph( fbc );
      else if ( ofs::magnetic_pole == 2 )
        Rotating_Dipole_LogSph( fbc );
    }
    { auto& fbc = fieldBC[UPPER_1];
      fbc.type = FieldBCType::DAMPING;
      fbc.damping_rate = ofs::damping_rate;
      fbc.indent = ofs::indent[UPPER_1];
    }
    { auto& fbc = fieldBC[LOWER_2];
      fbc.type = FieldBCType::COORDINATE;
      fbc.indent = ofs::indent[LOWER_2];
    }
    { auto& fbc = fieldBC[UPPER_2];
      fbc.type = FieldBCType::COORDINATE;
      fbc.indent = ofs::indent[UPPER_2];
    }
  }

  template < typename Real, int DGrid, typename RealJ >
  Updater<Real, DGrid, RealJ>::Updater( const mpi::CartComm& cart,
                                        const mani::Grid<Real,DGrid>& local_grid,
                                        apt::array< apt::pair<bool>, DGrid > is_at_boundary,
                                        int guard,
                                        Real (*omega_t) (Real),
                                        Real re_over_w_gyro_unitB )
    : _cart(cart),
      _4pi_times_re_over_w_gyro_unitB( 4.0 * std::acos(-1.0l) * re_over_w_gyro_unitB ) {
    static_assert( DGrid == 2 );
    for ( int i = 0; i < DGrid; ++i ) {
      fuparams.is_at_boundary[2*i] = is_at_boundary[i][LFT];
      fuparams.is_at_boundary[2*i + 1] = is_at_boundary[i][RGT];
      auto[ src, dest ] = cart.shift(i, 1);
      fuparams.neighbor_left[i] = src ? *(src) : NEIGHBOR_NULL;
      fuparams.neighbor_right[i] = dest ? *(dest) : NEIGHBOR_NULL;
    }
    for ( int i = DGrid; i < 3; ++i ) {
      fuparams.is_at_boundary[2*i] = false;
      fuparams.is_at_boundary[2*i + 1] = false;
      fuparams.neighbor_left[i] = NEIGHBOR_NULL;
      fuparams.neighbor_right[i] = NEIGHBOR_NULL;
    }

    // set fieldBC
    Set_FBC_RC_2DLogSph_Damp(fuparams.fieldBC, guard, omega_t);

    // grid
    { auto& grid = fuparams.grid;
      grid.dimension = DGrid;
      for ( int i = 0; i < DGrid; ++i ) {
        grid.dims[i] = local_grid[i].dim() + 2 * guard;
        grid.delta[i] = local_grid[i].delta();
        grid.lower[i] = local_grid[i].lower();

        grid.guard[i] = guard;
        for ( int j = 0; j < 2; ++j )
          grid.indent[2*i+j] = fuparams.fieldBC[2*i+j].indent;
      }
    }

    // NOTE Currently FieldUpdater only lives on primaries, which permits grid to be copied into Fields
    // set up E, B, J
    {
      const auto& grid = fuparams.grid;
      Efield.resize(grid);
      Bfield.resize(grid);
      current.resize(grid);
    }

    { // NOTE this is setting up E_bg and B_bg for damping and store in Efield, Bfield. IC for E and B are set up in operator(). This is to make resume work
      Efield.assign(0.0);
      Bfield.assign(0.0);
      const auto& grid = fuparams.grid;

      for (int j = 0; j < grid.dims[1]; j++) {
        Scalar theta = grid.pos(1, j, 0);
        Scalar theta_s = grid.pos(1, j, 1);
        for (int i = 0; i < grid.dims[0]; i++) {
          Scalar r = exp(grid.pos(0, i, 0));
          Scalar r_s = exp(grid.pos(0, i, 1));

          if ( ofs::magnetic_pole == 1 ) {
            Bfield(0, i, j ) = 1.0 / (r * r);
          } else if ( ofs::magnetic_pole == 2 ) {
            Bfield(0, i, j) = 2.0 * cos(theta_s) / (r * r * r);
            Bfield(1, i, j) = sin(theta) / (r_s * r_s * r_s);
          }
        }
      }

    }
    // create FieldUpdater
    fc.reset( new FieldCommunicator( cart, fuparams ) );
    fd.reset( new FiniteDiff(CoordType::LOG_SPHERICAL, fuparams.grid) );
    fu.reset( new FieldUpdater( fuparams, *fd, *fc, Efield, Bfield ) );

  }

  template < typename Real, int DGrid, typename RealJ >
  void Updater<Real,DGrid,RealJ>::operator() ( field_type& E,
                                               field_type& B,
                                               const J_type& Jmesh,
                                               Real dt, int timestep ) {
    static bool is_EB_initialized = false; // A hotfix to apply
    // NOTE
    // due to different stagger labeling systems, during conversion, only copy the bulk and send guards cells
    // For 0 <= i < dim_bulk,
    // 1) for unstaggered components (MIDWAY), indexing is field::Field[i] <-> VectorField[guard + i];
    // 2) for staggered components (INSITU), indexing is field::Field[i] <-> VectorField[guard + i - 1];

    { // convert from new to old
      const auto& grid = fuparams.grid;
      auto convert_from_new =
        [&]( auto& old_f, const auto& new_f ) {
          const auto& g = grid.guard;
          for ( int c = 0; c < 3; ++ c) {
            const auto& o = new_f[c].offset();
            for ( int j = 0; j < grid.reducedDim(1) + (o[1] == INSITU); ++j ) {
              for ( int i = 0; i < grid.reducedDim(0) + (o[0] == INSITU); ++i ) {
                old_f( c,
                       i + g[0] - (o[0] == INSITU),
                       j + g[1] - (o[1] == INSITU)
                       ) = new_f[c]({i,j});
              }}}
          // TODOL there is no need to copy guard cell values to current right? Because field solver doesn't use those values, instead they have boundary conditions
        };
      if ( !is_EB_initialized ) {
        convert_from_new( Efield, E );
        fc->SendGuardCells(Efield);
        convert_from_new( Bfield, B );
        fc->SendGuardCells(Bfield);
        is_EB_initialized = true;
      }
      convert_from_new( current, Jmesh );
      fc->SendGuardCells(current);

      // NOTE differences of new code
      // 1. The new code will evolve Maxwell's equations with 4\pi r_e / w_gyro.
      // 2. the new code passes in Jmesh, so conversion to real space current is needed
      for ( int c = 0; c < 3; ++c )
        for ( int i = 0; i < current.gridSize(); ++i )
          current.ptr(c)[i] *= _4pi_times_re_over_w_gyro_unitB;

      RestoreJToRealSpace(current, grid);
      fc->SendGuardCells(current);
    }

    {
      fu -> Update( Efield, Bfield, current, dt, timestep );
    }

    {// convert from old back to new
      auto convert_to_new =
        [&]( const auto& old_f, auto& new_f ) {
          const auto& grid = fuparams.grid;
          const auto& g = grid.guard;
          for ( int c = 0; c < 3; ++ c) {
            const auto& o = new_f[c].offset();
            for ( int j = 0; j < grid.reducedDim(1) + (o[1] == INSITU); ++j ) {
              for ( int i = 0; i < grid.reducedDim(0) + (o[0] == INSITU); ++i ) {
                new_f[c]({i,j}) = old_f( c,
                                         i + g[0] - (o[0] == INSITU),
                                         j + g[1] - (o[1] == INSITU)
                                         );
              }}}
        };
      convert_to_new( Efield, E );
      convert_to_new( Bfield, B );
      field::copy_sync_guard_cells( E, _cart );
      field::copy_sync_guard_cells( B, _cart );
    }
  }
}

namespace field::ofs {
  int magnetic_pole = 2;
  apt::array<int,4> indent {};
  double damping_rate = 10.0;
}

