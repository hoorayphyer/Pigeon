#include "CurrentDepositer.h"
#include "CoordSystem.h"
#include "mpi_wrapper.h"
#include "ShapeFunctions.h"
#include "Logger.h"
#include <chrono>
#include "InfoCollector.h"

#include "Domain.h"

namespace KahanJ {
  std::array<std::vector<Scalar>, 3> cJ;

  Scalar Sum( const Scalar& big_sum, const Scalar& small_input, Scalar& compensate ) {
    Scalar y = small_input - compensate;
    Scalar t = big_sum + y;
    compensate = ( t - big_sum ) - y;
    return t;
  }

  void Init( int size ) {
    for ( auto& c : cJ ) {
      c.resize(size);
      std::fill( c.begin(), c.end(), 0.0 );
    }
  }

}

inline double CalculateW_1DVersion( double sx0, double sx1 ){
    std::cout << "CalculateW_1DVersion not implemented yet!! " << std::endl;
    return 1.0;
}

inline double CalculateW_2DVersion( double sx0, double sx1, double sy0, double sy1 ){
    return ( 2 * sx1 * sy1 + sx0 * sy1 + sx1 * sy0 + 2 * sx0 * sy0 ) / 6.0;
}

inline double CalculateW_3DVersion( double sx0, double sx1,
                                    double sy0, double sy1,
                                    double sz0, double sz1 ){
    // return (   2 * ( sx1 * sy1 * sz1 - sx0 * sy1 * sz1 )
    //          + sx1 * sy0 * sz1 - sx0 * sy0 * sz1
    //          + sx1 * sy1 * sz0 - sx0 * sy1 * sz0
    //          + 2 * ( sx1 * sy0 * sz0 - sx0 * sy0 * sz0 )
    //         )/6.0;
  return (sx1 - sx0) * CalculateW_2DVersion(sy0, sy1, sz0, sz1);
}

struct DepositImpl {

  template< typename ShapeFunction >
  void operator() ( const ShapeFunction& shapefunction, VectorField<Scalar>& JField, double dt,
                    const Particles<QCarrier>& particles, const Grid& grid ) {
    auto radius = shapefunction.GetRadius();
    auto support = shapefunction.GetSupport();
    Scalar charge = particles.Attributes().charge;

    std::array<double,3> Wfactor = { -grid.delta[0]/dt, -grid.delta[1]/dt, -grid.delta[2]/dt };
    std::array<double,3> W = { 0.0, 0.0, 0.0 };

    for ( unsigned int index = 0; index < particles.Number(); index++) {
      if ( particles.IsEmpty( index ) ) continue;
      const auto& ptc = particles.PtcData()[index];
      if ( check_bit( ptc.flag, ParticleFlag::ignore_deposit ) ) continue;
      // use of initial position and cell. So this function
      // must be called before UpdatePos in Pusher.
      const auto& dx = ptc.dx;
      const auto& x0 = ptc.x;
      Vec3<int> cell = grid.getCell(ptc.cell);

      for( int k = - radius.z; k <= support.z - radius.z; k++){
        for( int j = - radius.y; j <=  support.y - radius.y; j++){
          for( int i = - radius.x; i <= support.x - radius.x; i++){
            int idx_ln = grid.getIdx( cell.x + i, cell.y + j, cell.z + k );
            double sx0 = shapefunction._sf1d( i + 0.5 - x0.x );
            double sy0 = shapefunction._sf1d( j + 0.5 - x0.y );

            double sx1 = shapefunction._sf1d( i + 0.5 - x0.x - dx.x);
            double sy1 = shapefunction._sf1d( j + 0.5 - x0.y - dx.y);


            //TODO: implement 1D.
            if( grid.dimension == 2 ){
              W[0] = CalculateW_3DVersion( sx0, sx1, sy0, sy1, 1.0, 1.0 );
              W[1] = CalculateW_3DVersion( sy0, sy1, 1.0, 1.0, sx0, sx1 );
              W[2] = CalculateW_2DVersion( sx0, sx1, sy0, sy1 );

              // FIXME fix the following
              // Calling deposition after pusher.calculateDisplacement implies
              // that p_tmp, which is at n+0.5 time step, is based with respect
              // to x^(n+1). However, the x used here is actually x^(n). One way
              // to improve this is obviously calling updatePos before deposition
              // and accordingly change expressions for calculating shapefunctions.

              // FIXME: But, where is J based? Does one really need rebasing momentum?
              const auto& p_tmp = ptc.p;
              Wfactor[2] = p_tmp.z / std::sqrt( 1.0 + p_tmp.dot(p_tmp) );
            } else if( grid.dimension == 3 ){
              double sz0 = shapefunction._sf1d( k + 0.5 - x0.z );
              double sz1 = shapefunction._sf1d( k + 0.5 - x0.z - dx.z);

              W[0] = CalculateW_3DVersion( sx0, sx1, sy0, sy1, sz0, sz1 );
              W[1] = CalculateW_3DVersion( sy0, sy1, sz0, sz1, sx0, sx1 );
              W[2] = CalculateW_3DVersion( sz0, sz1, sx0, sx1, sy0, sy1 );
            }

            // Note in the 2D case, JField(2,...) is really the value of J3, not
            // the difference. But we still store it in JField.

            for ( int i = 0; i < 3; ++i)
              JField.data(i) [idx_ln] = KahanJ::Sum( JField.data(i) [idx_ln], Wfactor[i] * W[i] * charge, KahanJ::cJ[i][idx_ln] );

          }
        }
      }
    }

    return;
  }

};

template < typename T >
void FoldBackMultiarrayAtAxis( MultiArray<T>& array, const Grid& grid, int axisDir,
                               int stag_axisDir, bool isupper, bool isadd = true ) {

  int trans[2] = { (axisDir + 1) % VECTOR_DIM , (axisDir + 2) % VECTOR_DIM };
  int increm[2] = { DirIncrement( trans[0], grid.extent() ) , DirIncrement( trans[1], grid.extent() ) };
  int increm_axis = DirIncrement( axisDir, grid.extent() );
  int istart = (isupper) ? grid.dims[axisDir]-2*grid.guard[axisDir] - isupper*stag_axisDir
    : grid.guard[axisDir]-stag_axisDir;
  for ( int k = 0; k < grid.dims[trans[0]]; k++ ) {
    for ( int j = 0; j < grid.dims[trans[1]]; j++ ) {
      int transIndex = k * increm[0] + j * increm[1];
      for ( int i = istart; i < istart + grid.guard[axisDir] + isupper*stag_axisDir; i++ ) {
        // find the mirror of i, whose value will be added to i
        int imirror = 2 * grid.guard[axisDir] - 1 - i - stag_axisDir
          + 2 * isupper * ( grid.dims[axisDir] - 2 * grid.guard[axisDir] );
        if ( isadd ) {
          array[ i * increm_axis + transIndex ] += array[ imirror * increm_axis + transIndex ];
        } else {
          array[ i * increm_axis + transIndex ] -= array[ imirror * increm_axis + transIndex ];
        }
        // erase folded-back contents
        // but when stag = 1, need to keep the values right on the boundaries,
        // which can be picked out by i == imirror
        if ( 0 == stag_axisDir || i != imirror )
          array[ imirror * increm_axis + transIndex ] = static_cast<T>( 0.0 );
      }
    }
  }

  return;
}

void HandleBoundaries( VectorField<Scalar>& JField, const Grid& grid, std::array<bool, 6> is_bdry, std::array<bool, 6> is_axis ) {
  for ( int i = 0; i < NUM_BOUNDARIES; ++i ) {
    if( ! is_bdry[i] ) continue;
    int axisDir = i / 2;
    bool isupper = i % 2;
    if ( grid.dimension == 2 && is_axis[i] ) {
      // Depending on when this function is called, JField can have different meanings,
      // which will affect the stagger and isadd when folding back.
      // Current implementation requires this function to be called before scanJ, which
      // means the JField here is really deltaRho. Therefore 0 is used in stagger_axis
      FoldBackMultiarrayAtAxis( JField.data(0), grid, axisDir, 0, isupper );
      FoldBackMultiarrayAtAxis( JField.data(1), grid, axisDir, 0, isupper );
      FoldBackMultiarrayAtAxis( JField.data(2), grid, axisDir, 0, isupper, false );
    }
  }
  return;
}

// the scheme works as follows. First sendadd transverse directions, then scan without
// tranverse guard cells, next sendadd the longitudinal direction, finally send guard
// cells of J, which is now done in RestoringJRhoInRealSpace.
void ScanJInCoordSpace(VectorField<Scalar>& JField, const Grid& grid, Domain& domain ){
  for ( int dir = 0; dir < VECTOR_DIM; dir++ ) {
    auto& Jcomp = JField.data(dir);
    int trans[2] = { ( dir + 1 ) % VECTOR_DIM, ( dir + 2 ) % VECTOR_DIM };

    // First: sendadd transverse directions
    //CLAIM: the case of sendAdding J3 in 2D is already properly accommodated.
    domain.SendAddCells( Jcomp, grid, trans[0] );
    domain.SendAddCells( Jcomp, grid, trans[1] );

    // Second: scan without transverse guard cells
    //Note that in the 2D case, no scan in J3 should be done. And this is automatically
    //fulfilled with the correct choice of grid.guard and grid.dims.
    int increm_trans[2] = { DirIncrement(trans[0], grid.extent()), DirIncrement(trans[1], grid.extent()) };
    int increm_long = DirIncrement( dir, grid.extent() );
    for ( int k = grid.guard[trans[0]]; k < grid.dims[trans[0]] - grid.guard[trans[0]]; k++ ) {
      for ( int j = grid.guard[trans[1]]; j < grid.dims[trans[1]] - grid.guard[trans[1]]; j++ ) {
        int transIndex = k * increm_trans[0] + j * increm_trans[1];
        Scalar compensation = 0.0;
        for ( int i = 1; i < grid.dims[dir]; i++ ) {
          // Jcomp[ i * increm_long + transIndex ] += Jcomp[ ( i - 1 ) * increm_long + transIndex ];
          Jcomp[ i * increm_long + transIndex ] = KahanJ::Sum( Jcomp[ ( i - 1 ) * increm_long + transIndex ], Jcomp[ i * increm_long + transIndex ], compensation );
        }
      }
    }

    // Third: sendadd longitudinal
    domain.SendAddCells( Jcomp, grid, dir, true );
  }
  return;
}

template< CoordType Coord >
void RestoreJToRealSpace( VectorField<Scalar>& JField, const Grid& grid ) {
  // static so that it can be used in lambda functions
  static typename CoordToScales<Coord>::type coord;

  const auto dim = grid.dimension;

  auto safe_divide =
    [] ( auto& n, auto d ) {
      if ( std::abs(d) > 1e-12 ) n /= d;
      else n = 0.0; };

  // define a function pointer.
  Scalar (*h_func) ( std::array<Scalar,3> q ) = nullptr;

  // FIXME TODO is the treatment also valid for 1D?
  //FIXME is it necessary to check the following for comp >= dim?
  // if( Coord == CoordType::LOG_SPHERICAL || Coord == CoordType::SPHERICAL
  //     || Coord == CoordType::LOG_SPHERICAL_EV ) {
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
          JField( comp, i, j, k ) /= h_func(q);
        }
      }
    }
  }

  return;
}


CurrentDepositer::CurrentDepositer( const Grid& grid ) { KahanJ::Init( grid.size() );}

template < CoordType Coord >
void CurrentDepositer::Deposit( VectorField<Scalar>& JField, double dt, const Particles<QCarrier>& ptcs, const AperParams& params, const MPICommunicator& comm, Domain& domain ) {
  const auto& grid = params.grid;
  auto weight_type = params.weightType;

  //FIXME TODO: zeroing should be done where?
  for ( auto& c : KahanJ::cJ ) {
    std::fill( c.begin(), c.end(), 0.0 );
  }
  // first deposit in coordinate space
  DepositImpl J_depositer;
  switch (grid.dimension) {
    // FIXME implement 1D.
    // case 1 :
    //   SelectWeight<1>( J_depositer, weight_type, JField, dt, ptcs, grid ); break;
  case 2 :
    SelectWeight<2>( J_depositer, weight_type, JField, dt, ptcs, grid ); break;
  case 3 :
    SelectWeight<3>( J_depositer, weight_type, JField, dt, ptcs, grid ); break;
  }

  HandleBoundaries( JField, grid, params.ens_specs.is_at_boundary, params.ens_specs.is_axis );

  comm.ensemble().barrier(); // for better timing
  auto t0 = high_resolution_clock::now();
  domain.EnsReduceFields( Jfield );
  auto t1 = high_resolution_clock::now();
  auto dur = duration_cast<clock_cycle>( t1 - t0 );
  InfoCollector::Instance().contents.t_ens_reduceJ = dur.count();

  if ( !comm.is_primary() ) return;

  ScanJInCoordSpace( JField, grid, domain );
  RestoreJRhoInRealSpace<Coord>( JField, grid );

  domain.SendGuardCells( JField );
}
