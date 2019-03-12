#ifndef  _DYNAMIC_VARIABLES_HPP_
#define  _DYNAMIC_VARIABLES_HPP_

#include "field/field.hpp"
#include "particle/array.hpp"

template< typename Real, std::size_t DGrid, std::size_t DPtc, typename state_t >
struct DynamicVars {
  field::Field<Real, 3, DGrid> E;
  field::Field<Real, 3, DGrid> B;
  field::Field<Real, 3, DGrid> J;

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
