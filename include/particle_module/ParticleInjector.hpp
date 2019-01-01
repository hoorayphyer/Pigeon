#ifndef PARTICLEINJECTOR_H
#define PARTICLEINJECTOR_H

#include "Dashboard.h"

enum class CoordType;
enum class WeightType;
class MPICommunicator;
struct AperData;
struct Rng;

class ParticleInjector {
private:
  DBPane_BoundaryCondition::PBC _pane;
  BoundaryPosition _bdry_pos;
  Scalar _Q_e;

  int NumToInject( int timestep, Scalar q1, Scalar q2, Scalar q3 ) const;

public:
  ParticleInjector( BoundaryPosition boundary_pos, Scalar Q_e, const DBPane_BoundaryCondition::PBC& pane );
  ~ParticleInjector();

  template < CoordType Coord >
  void InjectPairs(int timestep, Scalar dt, AperData& data, WeightType weight, const MPICommunicator& comm, Rng& rng );

};

#endif // PAIRMAKER_H
