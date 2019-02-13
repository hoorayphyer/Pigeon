#include "aperture.hpp"
#include "particle/updater.hpp"

template< typename Real, std::size_t DGrid, std::size_t DPtc, typename state_t >
Aperture<Real, DGrid, DPtc, state_t>::Aperture() {}

template< typename Real, std::size_t DGrid, std::size_t DPtc, typename state_t >
void Aperture<Real, DGrid, DPtc, state_t>::launch() {
  particle::Updater< Real, DGrid, DPtc, state_t,
                     knl::shape::Cloud_In_Cell,
                     particle::PairScheme::Disabled,
                     knl::coordsys_t::Cartesian,
                     particle::species::positron >
    update; // TODO
  for ( int i = _timestep_begin; i < _timestep_end; ++i ) {
    update( _dvars, _params, _grid, _comm );
  }
}
