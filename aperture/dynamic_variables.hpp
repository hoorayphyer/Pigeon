#ifndef  _DYNAMIC_VARIABLES_HPP_
#define  _DYNAMIC_VARIABLES_HPP_

#include "field/field.hpp"
#include "particle/map.hpp"
#include "particle/array.hpp"

template< typename Real, std::size_t DGrid, std::size_t DPtc, typename state_t >
struct DynamicVars {

  field::Field<Real, 3, DGrid> E;
  field::Field<Real, 3, DGrid> B;
  field::Field<Real, 3, DGrid> J;

  particle::map<particle::array<Real, DPtc, state_t>> particles;

  // ScalarField<Scalar> pairCreationEvents; // record the number of pair creation events in each cell.
  // PairCreationTracker pairCreationTracker;

};


#endif
