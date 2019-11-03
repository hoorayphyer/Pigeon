#include "field/old_field_solver/updater.hpp"

#include "field/field.hpp"
#include "field/sync.hpp"

#include "field/old_field_solver/FieldUpdater.h"
#include "field/old_field_solver/FieldCommunicator.h"
#include "field/old_field_solver/FiniteDiff.h"
#include "field/old_field_solver/CoordSystem.h"

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

  template < typename R, int DGrid, typename RJ >
  void OldSolve<R,DGrid,RJ>::Init( const mpi::CartComm& cart,
                                   const apt::Grid<R,DGrid>& local_grid) const {
    for ( int i = 0; i < DGrid; ++i ) {
      auto[ src, dest ] = cart.shift(i, 1);
      fuparams.neighbor_left[i] = src ? *(src) : NEIGHBOR_NULL;
      fuparams.neighbor_right[i] = dest ? *(dest) : NEIGHBOR_NULL;
    }
    for ( int i = DGrid; i < 3; ++i ) {
      fuparams.neighbor_left[i] = NEIGHBOR_NULL;
      fuparams.neighbor_right[i] = NEIGHBOR_NULL;
    }

    { // field BC, specialized for pulsar
      { auto& fbc = fuparams.fieldBC[LOWER_1];
        fbc.type = FieldBCType::ROTATING_CONDUCTOR;
        fbc.indent = _surface_indent;
        fbc.ft = _omega_t;
        if ( _magnetic_pole == 1 )
          Rotating_Monopole_LogSph( fbc );
        else if ( _magnetic_pole == 2 )
          Rotating_Dipole_LogSph( fbc );
      }
      { auto& fbc = fuparams.fieldBC[UPPER_1];
        fbc.type = FieldBCType::DAMPING;
        fbc.damping_rate = _damping_rate;
        fbc.indent = _damp_indent;
      }
      { auto& fbc = fuparams.fieldBC[LOWER_2];
        fbc.type = FieldBCType::COORDINATE;
        fbc.indent = _guard;
      }
      { auto& fbc = fuparams.fieldBC[UPPER_2];
        fbc.type = FieldBCType::COORDINATE;
        fbc.indent = _guard;
      }
    }

    // grid
    { auto& grid = fuparams.grid;
      grid.dimension = DGrid;
      for ( int i = 0; i < DGrid; ++i ) {
        grid.dims[i] = local_grid[i].dim() + 2 * _guard;
        grid.delta[i] = local_grid[i].delta();
        grid.lower[i] = local_grid[i].lower();

        grid.guard[i] = _guard;
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

          if ( _magnetic_pole == 1 ) {
            Bfield(0, i, j ) = 1.0 / (r * r);
          } else if ( _magnetic_pole == 2 ) {
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

  template < typename R, int DGrid, typename RJ >
  void OldSolve<R,DGrid,RJ>::operator() ( Field<R,3,DGrid>& E,
                                          Field<R,3,DGrid>& B,
                                          const Field<RJ,3,DGrid>& Jmesh,
                                          const apt::Grid<R,DGrid>& grid,
                                          const mpi::CartComm& cart,
                                          int timestep,
                                          R dt
                                          ) const {
    static bool is_initialized = false;
    if ( !is_initialized ) Init(cart, grid);

    // NOTE
    // due to different staggering systems, during conversion, only copy the bulk and send guards cells
    // For 0 <= i < dim_bulk,
    // 1) for unstaggered components (MIDWAY), indexing is field::Field[i] <-> VectorField[guard + i];
    // 2) for staggered components (INSITU), indexing is field::Field[i] <-> VectorField[guard + i - 1];

    { // convert from new to old
      const auto& fugrid = fuparams.grid;
      auto convert_from_new =
        [&]( auto& old_f, const auto& new_f ) {
          const auto& g = fugrid.guard;
          const auto beg_new = apt::range::begin(new_f.mesh().range());
          for ( int c = 0; c < 3; ++ c) {
            const auto& o = new_f[c].offset();
            for ( int j = 0; j < fugrid.reducedDim(1) + (o[1] == INSITU); ++j ) {
              for ( int i = 0; i < fugrid.reducedDim(0) + (o[0] == INSITU); ++i ) {
                old_f( c,
                       i + g[0] - (o[0] == INSITU),
                       j + g[1] - (o[1] == INSITU)
                       ) = new_f[c]({beg_new[0]+i,beg_new[1]+j});
              }}}
          // TODOL there is no need to copy guard cell values to current right? Because field solver doesn't use those values, instead they have boundary conditions
        };
      if ( !is_initialized ) {
        convert_from_new( Efield, E );
        fc->SendGuardCells(Efield);
        convert_from_new( Bfield, B );
        fc->SendGuardCells(Bfield);
        is_initialized = true;
      }
      convert_from_new( current, Jmesh );
      fc->SendGuardCells(current);

      // NOTE differences of new code
      // 1. The new code will evolve Maxwell's equations with 4\pi r_e / w_gyro.
      // 2. the new code passes in Jmesh, so conversion to real space current is needed
      for ( int c = 0; c < 3; ++c )
        for ( int i = 0; i < current.gridSize(); ++i )
          current.ptr(c)[i] *= _fourpi;

      RestoreJToRealSpace(current, fugrid);
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
          const auto beg_new = apt::range::begin(new_f.mesh().range());
          for ( int c = 0; c < 3; ++ c) {
            const auto& o = new_f[c].offset();
            for ( int j = 0; j < grid.reducedDim(1) + (o[1] == INSITU); ++j ) {
              for ( int i = 0; i < grid.reducedDim(0) + (o[0] == INSITU); ++i ) {
                new_f[c]({beg_new[0]+i,beg_new[1]+j}) = old_f( c,
                                                               i + g[0] - (o[0] == INSITU),
                                                               j + g[1] - (o[1] == INSITU)
                                                               );
              }}}
        };
      convert_to_new( Efield, E );
      convert_to_new( Bfield, B );
      field::copy_sync_guard_cells( E, cart );
      field::copy_sync_guard_cells( B, cart );
    }
  }
}
