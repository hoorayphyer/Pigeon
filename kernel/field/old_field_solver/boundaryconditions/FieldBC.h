#ifndef FIELDBC_H
#define FIELDBC_H

#include "Fields.h"
#include "Types.h"

class FieldBC {
 public:
  virtual void Apply(Scalar dt, VectorField<Scalar>& E, VectorField<Scalar>& B,
                     Scalar time) = 0;
  virtual ~FieldBC(){};
};

#include "FieldBC_coordinate.h"
#include "FieldBC_damping.h"
#include "FieldBC_rotating_conductor.h"

#endif  // FIELDBC_H
