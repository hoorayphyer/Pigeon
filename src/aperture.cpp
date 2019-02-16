#include "aperture.hpp"
#include "particle/updater.hpp"
#include "traits.hpp"

template< typename Real, std::size_t DGrid, std::size_t DPtc, typename state_t >
Aperture<Real, DGrid, DPtc, state_t>::Aperture() {}

template< typename Real, std::size_t DGrid, std::size_t DPtc, typename state_t >
void Aperture<Real, DGrid, DPtc, state_t>::launch() {
  particle::Updater< Real, DGrid, DPtc, state_t,
                     traits::shape,
                     traits::pair_produce_scheme,
                     traits::coordinate_system,
                     traits::posion_inj >
    update; // TODO
  for ( int i = _timestep_begin; i < _timestep_end; ++i ) {
    update( _dvars, _params, _grid, _ensemble );
  }
}

template class Aperture< traits::real_t, traits::DGrid, traits::DPtc, traits::ptc_state_t >;
