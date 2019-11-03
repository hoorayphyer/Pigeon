#include "FieldCommunicator.h"
#include "FiniteDiff.h"
#include "FieldUpdater.h"

#include "boundaryconditions/FieldBC.h"

#define UPPER_LIMIT 1.0e-10

static void CleanField (MultiArray<Scalar>& input, Scalar upper_lim) {
  for (int i = 0; i < input.size(); i++) {
    if (std::abs(input[i]) < upper_lim)
      input[i] = 0.0;
  }
}

template < typename FBCParams_t >
void SetFieldBC(std::unordered_map<BoundaryPosition, FieldBC*>& field_bdry_cond, const FUParams& params, const FBCParams_t& fieldBC, const VectorField<Scalar>& E_bg, const VectorField<Scalar>& B_bg ) {
  const auto isBoundary = params.is_at_boundary();
  for ( const auto& elm : fieldBC ) {
    auto b = elm.first;
    if ( !isBoundary[b] ) continue;

    const auto& fbc = elm.second;
    switch (fbc.type) {
    case FieldBCType::COORDINATE :
      field_bdry_cond[b] = new FieldBC_coordinate(b);
      break;
    case FieldBCType::DAMPING :
      field_bdry_cond[b] = new FieldBC_damping( b, params.grid, fieldBC.at(b), E_bg, B_bg );
      break;
    case FieldBCType::ROTATING_CONDUCTOR :
      field_bdry_cond[b] = new FieldBC_rotating_conductor( b, params.grid, fbc );
      break;
    } // end of switch
  } // end of for
}

FieldUpdater::FieldUpdater(const FUParams& params, FiniteDiff& fd, FieldCommunicator& fc, const VectorField<Scalar>& E_bg, const VectorField<Scalar>& B_bg)
  : _grid(params.grid), _is_at_boundary(params.is_at_boundary()),
    _fd(fd), _fc(fc) {

  const auto& grid = params.grid;

  SetFieldBC( _fieldBC, params, params.fieldBC, E_bg, B_bg );

  _BfieldOld = vector_field_type(grid);
  _BfieldDelta = vector_field_type(grid);
  _EfieldTmp = vector_field_type(grid);
  _EfieldDelta = vector_field_type(grid);

  _d_start = Index(grid.guard[0], grid.guard[1], grid.guard[2]);
  _d_ext = Extent(grid.reducedDim(0), grid.reducedDim(1), grid.reducedDim(2));

  // adjust _d_start and _d_ext. If at a true boundary, use indent instead of guard,
  // unless for the following boundary conditions: damping.

  for ( BoundaryPosition b = 0; b < NUM_BOUNDARIES; ++b ) {
    if ( !_is_at_boundary[b] ) continue;
    int dir = b / 2;
    bool islower = (b % 2 == 0);
    // check for damping
    if ( params.fieldBC.at(b).type == FieldBCType::DAMPING ) {
      continue;
    }
    else {
      _d_ext[dir] -= grid.indent[b] - grid.guard[dir];
      if (islower) {
        _d_start[dir] = grid.indent[b];
      }

    }
  }

  _fd.SetDefaultComputeDomain( _d_start, _d_ext);


  // Maybe the distinction between staggered and unstaggered fields?
  // Using the indent information to calculate the correct size of
  // domain of computation
//  for (int i = 0; i < NUM_BOUNDARIES; i++) {
//    if (!_domain -> isBoundary(i)) continue;
//    int dir = i / 2;
//    int indent = _domain -> params().indent[i];
//    if (indent > grid.guard[dir] && i % 2 == 0) {
//      _d_start[dir] = indent;
//      _d_ext[dir] -= (indent - grid.guard[dir]);
//    } else if (indent > grid.guard[dir]) {
//      _d_ext[dir] -= (indent - grid.guard[dir]);
//    }
//    // Is there anything wrong with the limits?
//    // if (pBC != nullptr && isConductor(pBC -> fieldBC((BoundaryPosition)i))) {
//    //   if (i % 2 == 0 && dir < grid.dim()) {
//    //     _d_start[dir] -= 1;
//    //     _d_ext[dir] += 1;
//    //   }
//    // }
//  }

}

FieldUpdater::~FieldUpdater() {
  for ( auto& fbc : _fieldBC )
    delete fbc.second;
}

void
FieldUpdater::ComputeBfieldUpdate( Scalar dt, VectorField<Scalar>& Efield, VectorField<Scalar>& Bfield, const VectorField<Scalar>& current) {
  // There is something wrong with the extent of the
  // laplacian. Maybe I neglected the distinction between staggered
  // and unstaggered fields

  _BfieldOld.copyFrom(Bfield);
  _fc.SendGuardCells(_BfieldOld);

  if (std::abs(_beta) > EPS) {
    _fd.ComputeLaplacian(Bfield, _BfieldDelta, FieldType::BTYPE, _is_at_boundary.data(), _fc, _d_start, _d_ext, true);

    _BfieldDelta.multiplyBy(_alpha * _beta * dt * dt);

    // for (int i = 0; i < VECTOR_DIM; i++)
    // CleanField(_BfieldDelta.data(i), UPPER_LIMIT);

    Bfield.addBy(_BfieldDelta);
  }

  _fc.SendGuardCells(Efield);


  // using _isBdry_EJ as an expediency. Should use _domain->isBoundary().
  bool isBdry_E[NUM_BOUNDARIES] = {};
  for ( int b = 0; b < NUM_BOUNDARIES; ++b ) {
    isBdry_E[b] = _is_at_boundary[b];
  }

  if ( _is_at_boundary[0] ) {
    isBdry_E[0] = false;
    _d_start[0] -= 1;
    _d_ext[0] += 1;
  }
  _fd.ComputeCurl(Efield, _BfieldDelta, FieldType::ETYPE, isBdry_E, _d_start, _d_ext);
  _BfieldDelta.multiplyBy(-1.0 * dt);
  Bfield.addBy(_BfieldDelta);

  if ( _is_at_boundary[0] ) {
    _d_start[0] += 1;
    _d_ext[0] -= 1;
    isBdry_E[0] = true;
  }

  _fd.ComputeCurl(current, _BfieldDelta, FieldType::ETYPE, isBdry_E, _d_start, _d_ext);
  _BfieldDelta.multiplyBy(_alpha * dt * dt);
  Bfield.addBy(_BfieldDelta);

  _fc.SendGuardCells(Bfield);

  _BfieldDelta.copyFrom(Bfield);

  Scalar factor = dt * dt * _alpha * _alpha;
  for (int j = 0; j < 5; ++j) {
    _BfieldDelta.multiplyBy(factor);
    _fd.ComputeLaplacian(_BfieldDelta, _BfieldDelta, FieldType::BTYPE, _is_at_boundary.data(), _fc, _d_start, _d_ext, true);
    // std::cout << "At iteration " << j << ": " << _BfieldDelta(0, 3, 257) << std::endl;
    // for (int i = 0; i < VECTOR_DIM; i++)
      // CleanField(_BfieldDelta.data(i), UPPER_LIMIT);
    Bfield.addBy(_BfieldDelta);
    _fc.SendGuardCells(_BfieldDelta);
  }

  _fc.SendGuardCells(Bfield);
}

void
FieldUpdater::ComputeBfieldUpdateExplicit(Scalar dt, VectorField<Scalar>& Efield, VectorField<Scalar>& Bfield, const VectorField<Scalar>& current) {
//  _domain->SendGuardCells(Efield);

  _fd.ComputeCurl(Efield, _BfieldDelta, FieldType::ETYPE, _is_at_boundary.data(), _d_start, _d_ext);
  _BfieldDelta.multiplyBy(-1.0 * dt);

  Bfield.addBy(_BfieldDelta);
  _fc.SendGuardCells(Bfield);
}

void
FieldUpdater::ComputeEfieldUpdate(Scalar dt, VectorField<Scalar>& Bfield, VectorField<Scalar>& Efield, const VectorField<Scalar>& current) {
  _EfieldDelta.assign(0.0);
//  _domain->SendGuardCells(Bfield);

  _fd.ComputeCurl(Bfield, _EfieldTmp, FieldType::BTYPE, _is_at_boundary.data(), _d_start, _d_ext);
  _EfieldTmp.multiplyBy(_alpha);
  _EfieldDelta.addBy(_EfieldTmp);

  _fd.ComputeCurl(_BfieldOld, _EfieldTmp, FieldType::BTYPE, _is_at_boundary.data(), _d_start, _d_ext);
  _EfieldTmp.multiplyBy(_beta);
  _EfieldDelta.addBy(_EfieldTmp);

  _EfieldDelta.subtractBy(current);
  _EfieldDelta.multiplyBy(dt);

  Efield.addBy(_EfieldDelta);

  for (int i = 0; i < VECTOR_DIM; i++)
    CleanField(Efield.data(i), UPPER_LIMIT);

  _fc.SendGuardCells(Efield);
}

void
FieldUpdater::ComputeEfieldUpdateExplicit(Scalar dt, VectorField<Scalar>& Bfield, VectorField<Scalar>& Efield, const VectorField<Scalar>& current) {
//  _domain->SendGuardCells(Bfield);

  _fd.ComputeCurl(Bfield, _EfieldDelta, FieldType::BTYPE, _is_at_boundary.data(), _d_start, _d_ext);
  _EfieldDelta.subtractBy(current);
  _EfieldDelta.multiplyBy(dt);
  Efield.addBy(_EfieldDelta);

  for (int i = 0; i < VECTOR_DIM; i++)
    CleanField(Efield.data(i), UPPER_LIMIT);

  _fc.SendGuardCells(Efield);
}

void
// see the line about elm.second -> Apply
FieldUpdater::Update(VectorField<Scalar>& Efield, VectorField<Scalar>& Bfield,
                     const VectorField<Scalar>& current, Scalar dt, int timestep) {

  ComputeBfieldUpdate(dt, Efield, Bfield, current);
  ComputeEfieldUpdate(dt, Bfield, Efield, current);
  // ComputeBfieldUpdateExplicit(dt, Efield, Bfield, current);
  // ComputeEfieldUpdateExplicit(dt, Bfield, Efield, current);

  // remove t_applyFieldBC from infocollector
  // need data rather than just E and B because damping will set particles flags. Improve this.
  for ( auto& elm : _fieldBC )
    elm.second->Apply( dt, Efield, Bfield, timestep * dt );

  // TODO TODO
  // auto& ssProxy = InfoCollector::Instance().ssProxy;
  // if ( ssProxy.saveThisSnapshot ) {
  //   CopyToSnapshot( ssProxy );
  // }

}
// TODO TODO
// // what is a way to remove this
// void FieldUpdater::CopyToSnapshot(SaveSnapshotProxy &proxy) const {
//   // the following BCs have data to be saved into snapshots: damping
//   for ( auto& elm : _fieldBC ) {
//     const auto bpos = elm.first;
//     auto& fbc = elm.second;
//     if ( _is_at_boundary[bpos] && fbc->bcType() == FieldBCType::DAMPING ) {
//       reinterpret_cast<FieldBC_damping*>(fbc) -> CopyToSnapshot(proxy);
//     }
//   }
// }
