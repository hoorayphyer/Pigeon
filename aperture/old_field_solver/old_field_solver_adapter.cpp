#include "old_field_solver_adapter.hpp"

#include "field/field.hpp"
#include "field/communication.hpp"

#include "FieldUpdater.h"
#include "FieldCommunicator.h"
#include "FiniteDiff.h"
#include "CoordSystem.h"

#include "parallel/mpi++.hpp"
#include <cmath>
#include <memory>

inline auto FT_SpinUp ( Scalar t_start, Scalar duration ) {
  return [t_start, duration](Scalar t) {
           Scalar t_diff = t - t_start;
           if ( t_diff < 0.0 ) return 0.0;
           else if (t_diff < duration ) return  t_diff / duration;
           else return 1.0;
         };
}

void Rotating_Monopole_LogSph( FBC& bdry, Scalar mu0, Scalar max_omega ) {
  bdry.B1 = [mu0] (Scalar r_log, Scalar , Scalar )
            {return mu0 / std::exp(2*r_log); };
  bdry.B2 = [] (Scalar, Scalar, Scalar) { return 0.0; };
  bdry.B3 = [] (Scalar, Scalar, Scalar){ return 0.0; };

  // find out E by E = -( Omega x r ) x B
  bdry.E1 = [=, B2=bdry.B2]( Scalar r_log, Scalar theta, Scalar phi ) {
              return max_omega * std::exp(r_log) * std::sin( theta ) * B2( r_log, theta, phi );};
  bdry.E2 = [=, B1=bdry.B1]( Scalar r_log, Scalar theta, Scalar phi ) {
              return - max_omega * std::exp(r_log) * std::sin( theta ) * B1( r_log, theta, phi );};
  bdry.E3 = []( Scalar, Scalar, Scalar ) {return 0.0;};
}

void Rotating_Dipole_LogSph( FBC& bdry, Scalar mu0, Scalar max_omega ) {
  bdry.B1 = [mu0] (Scalar x, Scalar theta, Scalar)
       { return mu0 * 2.0 * std::cos(theta) / std::exp(3*x); };
  bdry.B2 = [mu0] (Scalar x, Scalar theta, Scalar)
       { return mu0 * std::sin(theta) / std::exp(3*x); };
  bdry.B3 = [] (Scalar, Scalar, Scalar) { return 0.0; };

  // find out E by E = -( Omega x r ) x B
  bdry.E1 = [=, B2=bdry.B2]( Scalar r_log, Scalar theta, Scalar phi ) {
              return max_omega * std::exp(r_log) * std::sin( theta ) * B2( r_log, theta, phi );};
  bdry.E2 = [=, B1=bdry.B1]( Scalar r_log, Scalar theta, Scalar phi ) {
              return - max_omega * std::exp(r_log) * std::sin( theta ) * B1( r_log, theta, phi );};
  bdry.E3 = [=]( Scalar, Scalar, Scalar ) {return 0.0;};
}

namespace ofs {
  VectorField<Scalar> Efield;
  VectorField<Scalar> Bfield;
  VectorField<Scalar> current;

  FUParams fuparams;

  std::unique_ptr<FieldCommunicator> fc{};
  std::unique_ptr<FiniteDiff> fd{};
  std::unique_ptr<FieldUpdater> fu{};

  // Set field BCs for rotating conductor in 2D log spherical coordinates with a damping layer
  template < typename FBCs >
  void Set_FBC_RC_2DLogSph_Damp( FBCs& fieldBC, int guard ) {
    Scalar mu0 = 100.0;
    Scalar omega_max = 1.0 / 6.0;

    int magnetic_pole = 2; // 1 for mono-, 2 for di-
    int spinup_duration = 10.0;
    // field BC, specialized for pulsar
    { auto& fbc = fieldBC[LOWER_1];
      fbc.type = FieldBCType::ROTATING_CONDUCTOR;
      fbc.indent = 5;
      fbc.ft = FT_SpinUp(0.0, spinup_duration);
      if ( magnetic_pole == 1 )
        Rotating_Monopole_LogSph( fbc, mu0, omega_max );
      else if ( magnetic_pole == 2 )
        Rotating_Dipole_LogSph( fbc, mu0, omega_max );
    }
    { auto& fbc = fieldBC[UPPER_1];
      fbc.type = FieldBCType::DAMPING;
      fbc.damping_rate = 0.1;
      fbc.indent = 86;
    }
    { auto& fbc = fieldBC[LOWER_2];
      fbc.type = FieldBCType::COORDINATE;
      fbc.indent = guard;
    }
    { auto& fbc = fieldBC[UPPER_2];
      fbc.type = FieldBCType::COORDINATE;
      fbc.indent = guard;
    }
  }

  template < int DGrid >
  OldFieldUpdater<DGrid>::OldFieldUpdater( const mpi::CartComm& cart,
                                           const knl::Grid<double,DGrid,knl::grid1d::Clip>& local_grid,
                                           apt::array< apt::pair<bool>, DGrid > is_at_boundary,
                                           int guard) : _cart(cart) {
    for ( int i = 0; i < DGrid; ++i ) {
      fuparams.is_at_boundary[2*i] = is_at_boundary[i][LFT];
      fuparams.is_at_boundary[2*i + 1] = is_at_boundary[i][RGT];
      auto[ src, dest ] = cart.shift(i, 1);
      fuparams.neighbor_left[i] = src ? *(src) : NEIGHBOR_NULL;
      fuparams.neighbor_right[i] = dest ? *(dest) : NEIGHBOR_NULL;
    }
    for ( int i = DGrid; i < 3; ++i ) {
      fuparams.is_at_boundary[2*i] = true;
      fuparams.is_at_boundary[2*i + 1] = true;
      fuparams.neighbor_left[i] = NEIGHBOR_NULL;
      fuparams.neighbor_right[i] = NEIGHBOR_NULL;
    }

    // set fieldBC
    Set_FBC_RC_2DLogSph_Damp(fuparams.fieldBC, guard);

    // grid
    { auto& grid = fuparams.grid;
      grid.dimension = DGrid;
      for ( int i = 0; i < DGrid; ++i ) {
        grid.dims[i] = local_grid[i].dim() + 2 * guard;
        grid.delta[i] = local_grid[i].delta();
        grid.lower[i] = local_grid[i].lower(); // NOTE mesh containing this local_grid should have zero margin.

        grid.guard[i] = guard;
        for ( int j = 0; j < 2; ++j )
          grid.indent[2*i+j] = fuparams.fieldBC[2*i+j].indent;
      }
    }

    // NOTE Currently FieldUpdater only lives on primaries, which permits grid to be copied into Fields
    // set up E, B, J
    { auto& grid = fuparams.grid;
      Efield.resize(grid);
      Bfield.resize(grid);
      current.resize(grid);
    }

    // create FieldUpdater
    fc.reset( new FieldCommunicator( cart, fuparams ) );
    fd.reset( new FiniteDiff(CoordType::LOG_SPHERICAL, fuparams.grid) );
    fu.reset( new FieldUpdater( fuparams, *fd, *fc ) );
  }

  template < int DGrid >
  void OldFieldUpdater<DGrid>::operator() ( field_type& E,
                                            field_type& B,
                                            const field_type& J,
                                            double dt, int timestep ) {
    fuparams.dt = dt;
    // NOTE
    // due to different stagger labeling systems, during conversion, only copy the bulk and send guards cells
    // For 0 <= i < dim_bulk,
    // 1) for unstaggered components (MIDWAY), indexing is field::Field[i] <-> VectorField[guard + i];
    // 2) for staggered components (INSITU), indexing is field::Field[i] <-> VectorField[guard + i - 1];

    { // convert from new to old
      auto convert_from_new =
        [&]( auto& old_f, const auto& new_f ) {
          const auto& grid = fuparams.grid;
          const auto& g = grid.guard;
          for ( int c = 0; c < 3; ++ c) {
            const auto& o = new_f[c].offset();
            for ( int j = 0; j < grid.reducedDim(1); ++j ) {
              for ( int i = 0; i < grid.reducedDim(0); ++i ) {
                old_f( c,
                       i + g[0] - (o[0] == INSITU),
                       j + g[1] - (o[1] == INSITU)
                       ) = new_f[c]({i,j});
              }}}
          fc->SendGuardCells(old_f);
        };
      convert_from_new( Efield, E );
      convert_from_new( Bfield, B );
      convert_from_new( current, J );
      // NOTE
      // The new code will evolve Maxwell's equations with 4\pi. So `current` should first be multiplied by 4\pi
      for ( int c = 0; c < 3; ++c )
        for ( int i = 0; i < current.gridSize(); ++i )
          current.ptr(c)[i] *= (4 * CONST_PI);
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
            for ( int j = 0; j < grid.reducedDim(1); ++j ) {
              for ( int i = 0; i < grid.reducedDim(0); ++i ) {
                new_f[c]({i,j}) = old_f( c,
                                         i + g[0] - (o[0] == INSITU),
                                         j + g[1] - (o[1] == INSITU)
                                         );
                field::sync_guard_cells_from_bulk( new_f, _cart );
              }}}
        };
      convert_to_new( Efield, E );
      convert_to_new( Bfield, B );
    }
  }
}

namespace ofs {
  template struct OldFieldUpdater<>;
}
