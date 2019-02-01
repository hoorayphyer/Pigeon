#include "aperture.hpp"
#include "particle/updater.hpp"

template< typename Real_t, std::size_t DGrid, std::size_t DPtc >
Aperture::Aperture() {}

template< typename Real_t, std::size_t DGrid, std::size_t DPtc >
void Aperture::launch() {
  for ( int i = _timestep_i; i < _timestep_end; ++i ) {
    particle::update( _dvars, _params, _grid, _comm );
  }
}
