#ifndef _FIELDUPDATER_H_
#define _FIELDUPDATER_H_

#include <unordered_map>

#include "FUParams.h"
#include "Fields.h"

class FiniteDiff;
// struct SaveSnapshotProxy;
class FieldBC;
class FieldCommunicator;
////////////////////////////////////////////////////////////////////////////////

// TODO TODO put this part in the particle updater as a boundary condition
// // don't push particles that already entered damping layer
// const Grid& grid_full = data.Efield.grid();
// for ( auto& elm : data.particles ) {
//   auto& particles = elm.second;
//   // ignore neutral particles
//   if ( std::abs( particles.Attributes().charge ) < 1e-10 ) continue;

//   for ( int idx = 0; idx < particles.Number(); ++idx ) {
//     auto& ptc = particles.PtcData()[idx];
//     const int& C_dir = grid_full.getCell( ptc.cell )[dir];
//     // check if the particle lies inside the layer
//     if ( C_dir >= _shift[dir] && C_dir < _shift[dir] + _grid_damp.dims[dir] )
//     {
//       set_bit( ptc.flag, ParticleFlag::ignore_force );
//     }
//   }

// TODO TODO during restart, damping bg should be read in. There is the storing
// damping data into snapshot in FieldUpdater.cpp const auto& proxy =
// InfoCollector::Instance().ssProxy; if ( proxy.isReadDampingBg ) {
//   ReadFromSnapshot(proxy);
// } else {
// InitBackground( Efield_bg, Bfield_bg );
// }

// }
class FieldUpdater {
 public:
  typedef VectorField<Scalar> vector_field_type;
  typedef ScalarField<Scalar> scalar_field_type;

  FieldUpdater(const FUParams& params, FiniteDiff& fd, FieldCommunicator& fc,
               const VectorField<Scalar>& E_bg,
               const VectorField<Scalar>& B_bg);
  ~FieldUpdater();

  void Update(VectorField<Scalar>& Efield, VectorField<Scalar>& Bfield,
              const VectorField<Scalar>& current, Scalar dt, int timestep);

  // void CopyToSnapshot( SaveSnapshotProxy& proxy ) const;

 private:
  void ComputeEfieldUpdate(Scalar dt, VectorField<Scalar>& B,
                           VectorField<Scalar>& E,
                           const VectorField<Scalar>& current);
  void ComputeBfieldUpdate(Scalar dt, VectorField<Scalar>& E,
                           VectorField<Scalar>& B,
                           const VectorField<Scalar>& current);
  void ComputeEfieldUpdateExplicit(Scalar dt, VectorField<Scalar>& B,
                                   VectorField<Scalar>& E,
                                   const VectorField<Scalar>& current);
  void ComputeBfieldUpdateExplicit(Scalar dt, VectorField<Scalar>& E,
                                   VectorField<Scalar>& B,
                                   const VectorField<Scalar>& current);

  const Grid& _grid;
  std::array<bool, 6> _is_at_boundary;

  FiniteDiff& _fd;
  FieldCommunicator& _fc;

  Scalar _alpha = 1.0, _beta = 0.0;

  // contains fieldBC objects if the current ensemble is at physical boundary
  std::unordered_map<BoundaryPosition, FieldBC*> _fieldBC;

  Index _d_start;  // start of the computational domain on this node
  Extent _d_ext;   // extent of the computational domain on this node

  vector_field_type _BfieldOld;
  vector_field_type _BfieldDelta;
  vector_field_type _EfieldTmp;
  vector_field_type _EfieldDelta;

  //  // this is a hack for invoking ordinary diff along r for E and J
  //  transverse bool _isBdry_EJ[NUM_BOUNDARIES];

};  // ----- end of class FieldUpdater -----

#endif  // ----- #ifndef _FIELDUPDATER_H_  -----
