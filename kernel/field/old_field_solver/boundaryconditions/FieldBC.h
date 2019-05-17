#ifndef FIELDBC_H
#define FIELDBC_H

#include "Types.h"
#include "Fields.h"

class FieldBC {
public:
  virtual void Apply ( VectorField<Scalar>& E, VectorField<Scalar>& B, Scalar time ) = 0;
};

#include "FieldBC_coordinate.h"
#include "FieldBC_damping.h"
#include "FieldBC_rotating_conductor.h"

#endif // FIELDBC_H
