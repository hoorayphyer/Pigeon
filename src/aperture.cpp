#include "aperture.hpp"
#include "particle/updater.hpp"

template< typename Real, std::size_t DGrid, std::size_t DPtc >
Aperture<Real, DGrid, DPtc>::Aperture() {}

template< typename Real, std::size_t DGrid, std::size_t DPtc >
void Aperture<Real, DGrid, DPtc>::launch() {
  particle::Updater< Real, DGrid, DPtc,
                     knl::shape::Cloud_In_Cell,
                     particle::PairScheme::Disabled,
                     knl::coordsys_t::Cartesian,
                     particle::species::positron >
    update; // TODO
  for ( int i = _timestep_begin; i < _timestep_end; ++i ) {
    update( _dvars, _params, _grid, _comm );
  }
}
