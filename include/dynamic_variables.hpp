#ifndef  _DYNAMIC_VARIABLES_HPP_
#define  _DYNAMIC_VARIABLES_HPP_

#include "core/field.hpp"
#include "core/particle.hpp"

struct DynamicVars {
  // void init ( const DBPane_Particles& pane, PairProductionScheme pairScheme) {
  void init ( ) {
    // for ( auto[sp, num_max] : pane.capacity ) {
    //   particles[sp].resize(num_max);
    //   particles[sp].shrink_to_fit();
    // }

    // if (PairProductionScheme::DISABLED != pairScheme)
    //   pairCreationEvents.resize(grid);

  }

  Field<Real, 3, 3> E;
  Field<Real, 3, 3> B;
  Field<Real, 3, 3> j;

  species_map< std::vector<Particle> > particles;

  // ScalarField<Scalar> pairCreationEvents; // record the number of pair creation events in each cell.
  // PairCreationTracker pairCreationTracker;

};

#endif
