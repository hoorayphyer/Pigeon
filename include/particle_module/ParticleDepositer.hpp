#ifndef _PARTICLEDEPOSITER_H_
#define _PARTICLEDEPOSITER_H_

#include "Fields.h"

enum class CoordType;
class MPICommunicator;
struct Grid;
struct AperParams;

class Domain;

class ParticleDepositer {
  // FIXME TODO zeroing field should be done outside
  // FIXME TODO kahan summation not correct
  template < CoordType Coord, typename Ptc, typename Field, typename VarT >
  void Deposit( Field& field, const Particles<Ptc>& ptcs, const AperParams& params, const MPICommunicator& comm, Domain& domain, VarT (*f_get_var_from_ptc) (const Ptc&) );
};

#endif
