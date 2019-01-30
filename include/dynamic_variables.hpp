#ifndef  _DYNAMIC_VARIABLES_HPP_
#define  _DYNAMIC_VARIABLES_HPP_

#include "field/field.hpp"
#include "particle/particle.hpp"

template< typename T, std::size_t DGrid, std::size_t DPtc >
struct DynamicVars {
  Field<T, 3, DPtc> E;
  Field<T, 3, DPtc> B;
  Field<T, 3, DPtc> j;

  using particle::species;
  template < species sp >
  particle::vector<T, DPtc> particles;

  // ScalarField<Scalar> pairCreationEvents; // record the number of pair creation events in each cell.
  // PairCreationTracker pairCreationTracker;

};

namespace particle {
  constexpr auto fetch =
    []( species sp, auto& dynavars ) noexcept {
      switch ( sp ) {
      case species::electron : return dynavars.particles<species::electron>;
      case species::positron : return dynavars.particles<species::positron>;
      case species::ion : return dynavars.particles<species::ion>;
      case species::photon : return dynavars.particles<species::photon>;
      }
    };
}


#endif
