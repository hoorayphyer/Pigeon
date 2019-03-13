#include "old_field_solver_adapter.hpp"

#include "field/field.hpp"

#include "FieldUpdater.h"

namespace ofsa {
  VectorField<Scalar> Efield;
  VectorField<Scalar> Bfield;
  VectorField<Scalar> current;

  FUParams fuparams;
  Grid grid;

  OldFieldUpdater::OldFieldUpdater() {
    
  }

  OldFieldUpdater::~OldFieldUpdater() {
    
  }


  // TODO factor of 4\pi on J ?
  void OldFieldUpdater::operator( field_type& E,
                                  field_type& B,
                                  const field_type& J,
                                  field_type::element_t dt ) {
    // covert from new to old

    // do stuff here

    // covert from old back to new
  }
}
