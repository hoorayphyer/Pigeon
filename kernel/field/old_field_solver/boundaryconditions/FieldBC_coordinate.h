#ifndef FIELDBC_COORDINATE_H
#define FIELDBC_COORDINATE_H

#include "FieldBC.h"

class FieldBC_coordinate : public FieldBC {
private:
  const BoundaryPosition _bpos;

  template < typename T >
  void ReflectArrayAtAxis( MultiArray<T>& array, const Grid& grid, FieldType field_type, int comp, bool ispos = true ) const {

    int boundary = static_cast<int>(_bpos);
    int axisDir = boundary / 2;
    bool isupper = boundary % 2;
    int stag_axisDir = GetStagProperty( field_type, comp )[axisDir];
    int trans[2] = { (axisDir + 1) % VECTOR_DIM , (axisDir + 2) % VECTOR_DIM };
    int increm[2] = { DirIncrement( trans[0], grid.extent() ) , DirIncrement( trans[1], grid.extent() ) };
    int increm_axis = DirIncrement( axisDir, grid.extent() );
//    int istart = isupper ? grid.dims[axisDir] - 2 * grid.guard[axisDir] - isupper * stag_axisDir
//                           : grid.guard[axisDir] - stag_axisDir;
    int istart = isupper ? grid.dims[axisDir] - grid.indent[boundary] - stag_axisDir  : 0;
    int iend = isupper ? grid.dims[axisDir] : grid.indent[boundary];
    for ( int k = 0; k < grid.dims[trans[0]]; k++ ) {
      for ( int j = 0; j < grid.dims[trans[1]]; j++ ) {
        int transIndex = k * increm[0] + j * increm[1];
        for ( int i = istart; i < iend; i++ ) {
          // find the mirror of i
          int imirror = isupper ? 2 * grid.dims[axisDir] - 2 * grid.indent[boundary] - 1 - i - stag_axisDir
                                : 2 * grid.indent[boundary] - 1 - i - stag_axisDir;

          // if i = imirror and ispos = false, the value should be set to zero.
          if ( !ispos && i == imirror ) {
            array[ i * increm_axis + transIndex ] = 0.0;
          } else {
            array[ i * increm_axis + transIndex ] = array[ imirror * increm_axis + transIndex ];
            if (!ispos) array[ i * increm_axis + transIndex ] *= -1;
          }

        }
      }
    }

    return;
  }


public:
  FieldBC_coordinate( BoundaryPosition bpos ) : _bpos(bpos) {}

  virtual void Apply(Scalar dt, VectorField<Scalar>& Efield, VectorField<Scalar>& Bfield, Scalar time) override {
    const Grid& grid = Efield.grid();
    ReflectArrayAtAxis( Efield.data(0), grid, FieldType::ETYPE, 0 );
    ReflectArrayAtAxis( Efield.data(1), grid, FieldType::ETYPE, 1, false );
    ReflectArrayAtAxis( Efield.data(2), grid, FieldType::ETYPE, 2, false );

    ReflectArrayAtAxis( Bfield.data(0), grid, FieldType::BTYPE, 0 );
    ReflectArrayAtAxis( Bfield.data(1), grid, FieldType::BTYPE, 1, false );
    ReflectArrayAtAxis( Bfield.data(2), grid, FieldType::BTYPE, 2, false );
  }

};

#endif // FIELDBC_COORDINATE_H
