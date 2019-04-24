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

#include "gen.hpp"

inline auto FT_SpinUp (Scalar t) noexcept {
  return std::min(t / pic::spinup_duration, 1.0);
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

void RestoreJToRealSpace( VectorField<Scalar>& JField, const Grid& grid ) {
  // static so that it can be used in lambda functions
  static typename CoordToScales<CoordType::LOG_SPHERICAL>::type coord;

  const auto dim = grid.dimension;

  // define a function pointer.
  Scalar (*h_func) ( std::array<Scalar,3> q ) = nullptr;

  for ( int comp = 0; comp < dim; ++comp ) {
    Index stagJ;
    // NOTE only non-capture lambda can be assigned to function pointer
    if ( comp < dim ) {
      stagJ = GetStagProperty(FieldType::ETYPE, comp);
      switch(comp) {
      case 0 : h_func =
          [] ( std::array<Scalar,3> q ) {
            return coord.h2( q[0], q[1], q[2] ) * coord.h3( q[0], q[1], q[2] ); }; break;
      case 1 : h_func =
          [] ( std::array<Scalar,3> q ) {
            return coord.h3( q[0], q[1], q[2] ) * coord.h1( q[0], q[1], q[2] ); }; break;
      case 2 : h_func =
          [] ( std::array<Scalar,3> q ) {
            return coord.h1( q[0], q[1], q[2] ) * coord.h2( q[0], q[1], q[2] ); }; break;
      }
    } else {
      stagJ = Index(0, 0, 0);
      h_func =
        [] ( std::array<Scalar,3> q ) {
          return coord.h1( q[0], q[1], q[2] ) * coord.h2( q[0], q[1], q[2] )
            * coord.h3( q[0], q[1], q[2] ); };
    }

    std::array<Scalar,3> q = { 0, 0, 0 };
    for ( int k = grid.guard[2]; k < grid.dims[2] - grid.guard[2]; ++k ) {
      q[2] = grid.pos( 2, k, stagJ[2] );
      for ( int j = grid.guard[1]; j < grid.dims[1] - grid.guard[1]; ++j ) {
        q[1] = grid.pos( 1, j, stagJ[1] );
        for( int i = grid.guard[0]; i < grid.dims[0] - grid.guard[0]; ++i ) {
          q[0] = grid.pos( 0, i, stagJ[0]);

          // FIXME is it OK to leave out check on h_func(q) being nearly zero?
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
    // field BC, specialized for pulsar
    { auto& fbc = fieldBC[LOWER_1];
      fbc.type = FieldBCType::ROTATING_CONDUCTOR;
      fbc.indent = pic::ofs::indent[LOWER_1];
      fbc.ft = FT_SpinUp;
      if ( pic::ofs::magnetic_pole == 1 )
        Rotating_Monopole_LogSph( fbc, pic::mu0, pic::omega_max );
      else if ( pic::ofs::magnetic_pole == 2 )
        Rotating_Dipole_LogSph( fbc, pic::mu0, pic::omega_max );
    }
    { auto& fbc = fieldBC[UPPER_1];
      fbc.type = FieldBCType::DAMPING;
      fbc.damping_rate = pic::ofs::damping_rate;
      fbc.indent = pic::ofs::indent[UPPER_1];
    }
    { auto& fbc = fieldBC[LOWER_2];
      fbc.type = FieldBCType::COORDINATE;
      fbc.indent = pic::ofs::indent[LOWER_2];
    }
    { auto& fbc = fieldBC[UPPER_2];
      fbc.type = FieldBCType::COORDINATE;
      fbc.indent = pic::ofs::indent[UPPER_2];
    }
  }

  template < int DGrid >
  OldFieldUpdater<DGrid>::OldFieldUpdater( double unit_e,
                                           const mpi::CartComm& cart,
                                           const knl::Grid<double,DGrid>& local_grid,
                                           apt::array< apt::pair<bool>, DGrid > is_at_boundary,
                                           int guard ) : _cart(cart), _unit_e(unit_e) {
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
                                            const J_type& Jmesh,
                                            double dt, int timestep ) {
    fuparams.dt = dt;
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
      convert_from_new( current, Jmesh );

      // NOTE differences of new code
      // 1. The new code will evolve Maxwell's equations with 4\pi. So `current` should first be multiplied by 4\pi
      // 2. the new code passes in Jmesh, so conversion to real space current is needed
      // 3. Jmesh needs to be multiplied by unit_e
      for ( int c = 0; c < 3; ++c )
        for ( int i = 0; i < current.gridSize(); ++i )
          current.ptr(c)[i] *= (4 * CONST_PI * _unit_e);

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
