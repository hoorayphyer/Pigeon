#include "ParticleInjector.h"
#include "CoordSystem.h"
#include "Logger.h"
#include "AperData.h"
#include "mpi_wrapper.h"
#include <ctime>
#include <random>
#include <algorithm>
#include <stdexcept>
#include "cmath"
#include "ShapeFunctions.h"
#include "InfoCollector.h"
#include "Rng.h"

ParticleInjector::ParticleInjector( BoundaryPosition boundary_pos, Scalar Q_e, const DBPane_BoundaryCondition::PBC& pane )
  : _bdry_pos(boundary_pos), _Q_e(Q_e), _pane(pane) {}

ParticleInjector::~ParticleInjector() {}

int ParticleInjector::NumToInject(int timestep, Scalar q1, Scalar q2, Scalar q3) const {
  // the idea is that (timestep + 1) * N_inj * profile represents the accummulated number of injected pairs through the specified timestep. NOTE the timestep is shifted by one to reflect the actual times the inject is called, including this time.
  Scalar inj_num_base = _pane.N_inj * _pane.x_profile( q1, q2, q3 );
  return static_cast<int>( (timestep + 1) * inj_num_base ) - static_cast<int>( timestep * inj_num_base );
}

// the number needed to be injected to sustain the current
inline Scalar N_j( Scalar jr, Scalar Q_e, Scalar theta, Scalar dt, Scalar dr, Scalar multiple ) {
  return multiple * jr * std::sin(theta) * dt / (Q_e * dr);
}

class MovingAverage {
private:
  std::vector<Scalar> _average;

public:
  MovingAverage( int num_averages ) {
    _average.resize(num_averages);
    std::fill( _average.begin(), _average.end(), 0.0 );
  }

  Scalar Average( Scalar f_current, int index, Scalar beta ) {
    // effectively the average is approximately over 1 / ( 1 - beta ) data points
    _average[index] = beta * _average[index] + ( 1 - beta ) * f_current;
    return _average[index];
  }
};

std::unique_ptr<MovingAverage> mov_av = nullptr;

inline void local_EB( Vec3<Scalar>& EVec, Vec3<Scalar>& BVec, const AperData& data, const Index& vecIdx, WeightType weight ) {
  static FieldInterpolater interpolator;
  const auto& EField = data.Efield;
  const auto& BField = data.Bfield;
  const int dim = EField.grid().dimension;
  Vec3<POS_TYPE> x_interp(0.5, 0.5, 0.5); // note shapefunction has taken care of dim, so 0.5 beyond _dim will be ignored automatically.
  if ( 2 == dim ) {
    SelectWeight<2>( interpolator, weight, EVec, EField, FieldType::ETYPE, vecIdx, x_interp );
    SelectWeight<2>( interpolator, weight, BVec, BField, FieldType::BTYPE, vecIdx, x_interp );
  } else if ( 3 == dim ) {
    SelectWeight<3>( interpolator, weight, EVec, EField, FieldType::ETYPE, vecIdx, x_interp );
    SelectWeight<3>( interpolator, weight, BVec, BField, FieldType::BTYPE, vecIdx, x_interp );
  }
}

template < CoordType Coord >
void ParticleInjector::InjectPairs( int timestep, Scalar dt, AperData& data, WeightType weight, const MPICommunicator& comm, Rng& rng ) {
  typename CoordToScales<Coord>::type coord;
  Scalar time = timestep * dt;

  const auto& ensemble = comm.ensemble();
  int ens_rank_inj = timestep % ensemble.size();

  if ( ensemble.rank() != ens_rank_inj ) return;

  int dir = _bdry_pos / 2;
  const auto& J = data.Jfield;
  const Grid& grid = J.grid();
  const int dim = grid.dimension;
  int trans[2] = { ( dir + 1 ) % VECTOR_DIM, ( dir + 2 ) % VECTOR_DIM };

  Index vecIdx;
  vecIdx[dir] = grid.indent[_bdry_pos] - 1; // inject right below the indent

  auto& ptcs_minus = data.particles.at(_pane.q_minus);
  auto& ptcs_plus = data.particles.at(_pane.q_plus);
  // auto& electron_tracker = electrons.Tracker();
  // auto& pos_charge_tracker = pos_charges.Tracker();

  if ( !mov_av ) {
    mov_av.reset( new MovingAverage( grid.reducedDim(trans[0]) * grid.reducedDim(trans[1]) ) );
    // set the mov_av to initial values of Jr, which is done by using a beta =
    // 0. Since we don't store mov_av into snapshot, this is to avoid a huge
    // plumit in averaged Jr when the run is resumed.
    for ( int k = grid.guard[trans[0]]; k < grid.dims[trans[0]] - grid.guard[trans[0]]; ++k ) {
      for ( int j = grid.guard[trans[1]]; j < grid.dims[trans[1]] - grid.guard[trans[1]]; ++j ) {
        vecIdx[trans[1]] = j;
        vecIdx[trans[0]] = k;
        mov_av -> Average( J(0,vecIdx[0], vecIdx[1], vecIdx[2]), vecIdx[trans[0]] + vecIdx[trans[1]] * grid.reducedDim(trans[0]), 0.0 );
      }
    }
  }

  // Here it is assumed that profile only covers the transverse directions.
  for ( int k = grid.guard[trans[0]]; k < grid.dims[trans[0]] - grid.guard[trans[0]]; ++k ) {
    for ( int j = grid.guard[trans[1]]; j < grid.dims[trans[1]] - grid.guard[trans[1]]; ++j ) {
      vecIdx[trans[1]] = j;
      vecIdx[trans[0]] = k;

      auto q1 = grid.pos( 0, vecIdx[0], 0 );
      // FIXME: the conditional may not be necessary because coordinates
      // corresponding to symmetry can be traced normally
      auto q2 = dim >= 2 ? grid.pos( 1, vecIdx[1], 0 ) : grid.lower[1];
      auto q3 = dim >= 3 ? grid.pos( 2, vecIdx[2], 0 ) : grid.lower[2];

      // To find rho_bound, find the longest edge of the physical volume element
      auto lambda_p = coord.h1( q1, q2, q3 ) * grid.delta[0];
      if ( dim >= 2 ) {
        auto tmp =  coord.h2( q1, q2, q3 ) * grid.delta[1];
        lambda_p = std::max( lambda_p, tmp );
      }
      if ( dim >= 3 ) {
        auto tmp =  coord.h3( q1, q2, q3 ) * grid.delta[2];
        lambda_p = std::max( lambda_p, tmp );
      }
      auto rho_bound = 0.8 * 4 * CONST_PI * CONST_PI / ( lambda_p * lambda_p );
      auto cell_size = coord.h1(q1, q2, q3) * coord.h2(q1, q2, q3) * coord.h3(q1, q2, q3);
      auto n_bound = rho_bound * cell_size / _Q_e;
      // below data.ptcNumbers are used which consist of local particles on the process. Assuming a reasonably balanced load, n_bound should be divided by the ensemble size
      n_bound /= ensemble.size();

      // calculate current value of n
      Scalar n_plus = 0.0;
      if ( data.particles.find(ParticleType::ION) != data.particles.end() ) {
        n_plus += data.ptcNumbers.at(ParticleType::ION) (vecIdx[0], vecIdx[1], vecIdx[2]);
      }
      if ( data.particles.find(ParticleType::POSITRON) != data.particles.end() ) {
        n_plus += data.ptcNumbers.at(ParticleType::POSITRON) (vecIdx[0], vecIdx[1], vecIdx[2]);
      }
      Scalar n_minus = data.ptcNumbers.at(ParticleType::ELECTRON) (vecIdx[0], vecIdx[1], vecIdx[2]);

      Vec3<Scalar> EVec, n_B;
      local_EB(EVec, n_B, data, vecIdx, weight);
      n_B /= std::sqrt( n_B.dot( n_B ) ); // TODO: add check on nonzeroness of n_B?
      n_B *= ( 0.0 <= n_B.x ) - ( n_B.x < 0.0 ); // make n_B always point away from the star
      Scalar n_current = ( EVec.dot(n_B) > 0 ) ? n_plus : n_minus;

      if ( n_bound > n_current ) {
        auto Jr = mov_av -> Average( J(0,vecIdx[0], vecIdx[1], vecIdx[2]), vecIdx[trans[0]] + vecIdx[trans[1]] * grid.reducedDim(trans[0]), std::min( 1.0, 1.0 - dt / _pane.t_moving_average ) );
        auto num_inj = N_j( Jr, _Q_e, q2, dt, grid.delta[0], _pane.J_reg_x );
        num_inj = std::max( num_inj, static_cast<Scalar>(NumToInject( timestep, q1, q2, q3)) );
        int num = static_cast<int>( std::min( num_inj, n_bound - n_current ) );
        InfoCollector::Instance().contents.num_injectedPairs += num * 2;

        int cell = grid.getIdx( vecIdx[0], vecIdx[1], vecIdx[2] );
        for ( int n = 0; n < num; ++n ) {
          Vec3<POS_TYPE> x;
          x[dir] = 1.0 - _pane.layer * rng.uniform();
          if ( trans[0] < dim )
            x[ trans[0] ] = rng.uniform();
          if ( trans[1] < dim )
            x[ trans[1] ] = rng.uniform();
          // sample momentum
          auto pos_full = grid.pos_particle( cell, x );
          auto p = _pane.p_profile( time, pos_full.x, pos_full.y, pos_full.z  );
          // add a random component along B field
          p += n_B * ( _pane.sigma_p.x * rng.gaussian() );
          QCarrier ptc;
          ptc.x = x, ptc.p = p, ptc.cell = cell;
          ptcs_minus.Append( ptc );
          ptcs_plus.Append( ptc );
          // // check tracking
          // if ( electron_tracker.IsDoTrack( timestep, rng.uniform() ) )
          //   electron_tracker.Track( electrons.PtcData()[electrons.Number() - 1] );
          // if ( pos_charge_tracker.IsDoTrack( timestep, rng.uniform() ) )
          //   pos_charge_tracker.Track( pos_charges.PtcData()[pos_charges.Number() - 1] );
        }
      }
    }
  }

}

// INSTANTIATE_FUNCTION_WITH_COORDTYPES(void, ParticleInjector::InjectPairs, ( int timestep, Scalar dt, AperData& data, WeightType weight, const MPICommunicator& comm, Rng& rng ) );
