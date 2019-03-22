#ifndef FIELDBC_H
#define FIELDBC_H

#include "Types.h"
#include "Fields.h"

class FieldBC {
public:
  virtual void Apply ( VectorField<Scalar>& E, VectorField<Scalar>& B, Scalar time ) = 0;
  virtual FieldBCType bcType() const = 0;
};

#include "boundaryconditions/FieldBC_coordinate.h"
#include "boundaryconditions/FieldBC_damping.h"
#include "boundaryconditions/FieldBC_rotating_conductor.h"

#endif // FIELDBC_H
