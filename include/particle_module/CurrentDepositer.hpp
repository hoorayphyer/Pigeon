#ifndef CURRENTDEPOSITER_H
#define CURRENTDEPOSITER_H

#include "Fields.h"
#include "Particles.h"

enum class CoordType;
class MPICommunicator;
struct Grid;
struct AperParams;

class Domain;

class CurrentDepositer{
public:
  CurrentDepositer( const Grid& grid );
  ~CurrentDepositer() {}

  // FIXME TODO zeroing field should be done outside
  // FIXME TODO kahan summation not correct
  template < CoordType Coord >
  void Deposit( VectorField<Scalar>& JField, double dt, const Particles<QCarrier>& ptcs, const AperParams& params, const MPICommunicator& comm, Domain& domain );

};

#endif
