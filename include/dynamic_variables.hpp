#ifndef  _DYNAMIC_VARIABLES_HPP_
#define  _DYNAMIC_VARIABLES_HPP_

#include "field/field.hpp"
#include "particle/array.hpp"

template< typename Real, int DGrid, int DPtc, typename state_t >
struct DynamicVars {
  field::Field<Real, 3, DPtc> E;
  field::Field<Real, 3, DPtc> B;
  field::Field<Real, 3, DPtc> J;

  particle::array<Real, DPtc, state_t> electrons;
  particle::array<Real, DPtc, state_t> positrons;
  particle::array<Real, DPtc, state_t> ions;
  particle::array<Real, DPtc, state_t> photons;

  constexpr auto& operator[] ( particle::species sp ) noexcept {
    using particle::species;
    switch ( sp ) {
    case species::electron : return electrons;
    case species::positron : return positrons;
    case species::ion : return ions;
    case species::photon : return photons;
    }
  }

  constexpr const auto& operator[] ( particle::species sp ) const noexcept {
    using particle::species;
    switch ( sp ) {
    case species::electron : return electrons;
    case species::positron : return positrons;
    case species::ion : return ions;
    case species::photon : return photons;
    }
  }

  // ScalarField<Scalar> pairCreationEvents; // record the number of pair creation events in each cell.
  // PairCreationTracker pairCreationTracker;

};


#endif
