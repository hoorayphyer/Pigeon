#include "ParticleDepositer.h"
#include "CoordSystem.h"
#include "mpi_wrapper.h"
#include "ShapeFunctions.h"
#include "Logger.h"
#include <chrono>
#include "InfoCollector.h"

#include "Domain.h"

struct Kahan {
  std::vector<Scalar> c; // compensations

  Scalar Sum( const Scalar& big_sum, const Scalar& small_input, Scalar& compensate ) {
    Scalar y = small_input - compensate;
    Scalar t = big_sum + y;
    compensate = ( t - big_sum ) - y;
    return t;
  }

  Kahan( int size ) {
    c.reserve(size);
    c.resize(size);
    std::fill( c.begin(), c.end(), 0.0 );
  }

};

struct DepositImpl {
  template < typename Field >
  struct VarT {};

  template <typename T>
  struct VarT<ScalarField<T>> { using type = T;};

  template <typename T>
  struct VarT<VectorField<T>> { using type = Vec3<T>;};

  template< typename Field, typename ShapeFunction, typename Ptc >
  void operator() ( const ShapeFunction& shapefunction, Field& field, const Particles<Ptc>& particles, const Grid& grid, typename VarT<Field>::type (*f) (const Ptc&) ) {
    // FIXME might need same kahan across species deposit? Or not, because it will only miss three compensates if there are three species.
    constexpr bool is_sf = std::is_same<Field, ScalarField<Scalar> >:: value;
    std::array< Kahan, is_sf ? 1 : 3 > kh = {grid.size()}; // FIXME check syntax. Initialize all Kahan with grid.size()

    auto radius = shapefunction.GetRadius();
    auto support = shapefunction.GetSupport();
    for ( unsigned int index = 0; index < particles.Number(); index++) {
      if ( particles.IsEmpty( index ) ) continue;
      const auto& ptc = particles.PtcData()[index];
      if ( check_bit( ptc.flag, ParticleFlag::ignore_deposit ) ) continue;
      const auto& dx = ptc.dx;
      const auto& x0 = ptc.x;
      const auto cell = grid.getCell(ptc.cell);

      auto var_deposit = f(ptc);
      for( int k = - radius.z; k <= support.z - radius.z; k++){
        for( int j = - radius.y; j <=  support.y - radius.y; j++){
          for( int i = - radius.x; i <= support.x - radius.x; i++){
            double weight = shapefunction._sf1d( i + 0.5 - x0.x - dx.x);
            if ( grid.dimension > 1 )
              weight *= shapefunction._sf1d( j + 0.5 - x0.y - dx.y);
            if( grid.dimension > 2 )
              weight *= shapefunction._sf1d( k + 0.5 - x0.z - dx.z);

            int idx_ln = grid.getIdx( cell.x + i, cell.y + j, cell.z + k );
            if constexpr
              ( is_sf )
                field.data()[idx_ln] = kh[0]( field[idx_ln], weight * var_deposit, kh[0]::c[idx_ln] );
            else  { // deal with vectorfield
              for (int i = 0; i < 3; ++i )
                field.data(i)[idx_ln] = kh[i]( field[idx_ln], weight * var_deposit[i], kh[i]::c[idx_ln] );
              }

          }
        }
      }
    }

    return;
  }

};

// FIXME this function also exists in currentDepositer
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

void HandleBoundaries( MultiArray<Scalar>& field, const Grid& grid, std::array<bool, 6> is_bdry, std::array<bool, 6> is_axis ) {
  for ( int i = 0; i < NUM_BOUNDARIES; ++i ) {
    if( ! is_bdry[i] ) continue;
    int axisDir = i / 2;
    bool isupper = i % 2;
    if ( grid.dimension == 2 && is_axis[i] ) {
      FoldBackMultiarrayAtAxis( field, grid, axisDir, 0, isupper );
    }
  }
  return;
}

template < CoordType Coord, typename Ptc, typename Field, typename VarT >
void ParticleDepositer::Deposit( Field& field, const Particles<Ptc>& ptcs, const AperParams& params, const MPICommunicator& comm, Domain& domain, VarT (*f) (const Ptc&) ) {
  constexpr bool is_sf = std::is_same<Field, ScalarField<Scalar> >:: value;
  const auto& grid = params.grid;
  auto weight_type = params.weightType;

  // first deposit in coordinate space
  DepositImpl particle_depositer;

  switch (grid.dimension) {
  case 1 :
    // FIXME check 1D.
    SelectWeight<1>( particle_depositer, weight_type, field, ptcs, grid, f ); break;
  case 2 :
    SelectWeight<2>( particle_depositer, weight_type, field, ptcs, grid, f ); break;
  case 3 :
    SelectWeight<3>( particle_depositer, weight_type, field, ptcs, grid, f ); break;
  }

  if constexpr
    ( is_sf )
      HandleBoundaries( field.data(), grid, params.ens_specs.is_at_boundary, params.ens_specs.is_axis );
  else {
    for ( int n = 0; n < 3; ++n )
      HandleBoundaries( field.data(n), grid, params.ens_specs.is_at_boundary, params.ens_specs.is_axis );
  }

  domain.EnsReduceFields( field );

  if ( !comm.is_primary() ) return;

  for ( int i = 0; i < grid.dimension; i++ )
    domain.SendAddCells( field, grid, i );

  // restore field to real space

  auto safe_div = [] ( auto& num, auto& den ) {
                    if ( std::abs(den) > 1e-12 ) {
                      num /= den;
                    } else { // if hhh = 0
                      num = 0.0;
                    }
                  };

  typename CoordToScales<Coord>::type coord;
  for ( int k = grid.guard[2]; k < grid.dims[2] - grid.guard[2]; ++k ) {
    auto q3 = grid.pos( 2, k, 0 );
    for( int j = grid.guard[1]; j < grid.dims[1] - grid.guard[1]; ++j ) {
      auto q2 = grid.pos( 1, j, 0 );
      for( int i = grid.guard[0]; i < grid.dims[0] - grid.guard[0]; ++i ) {
        auto q1 = grid.pos( 0, i, 0 );

        auto hhh = coord.h1(q1,q2,q3) * coord.h2(q1,q2,q3) * coord.h3(q1,q2,q3) * grid.dV;
        if constexpr
          ( is_sf ) {
            safe_div( field ( i, j, k ), hhh );
          } else {
          for ( int n = 0; n < 3; ++n )
            safe_div( field(n, i, j, k ), hhh );
        }
      }
    }
  }

  // send guard cells
  domain.SendGuardCells( field );
}
