#ifndef _PAIR_PRODUCER_HPP_
#define _PAIR_PRODUCER_HPP_

#include "types.hpp"
class Rng;

namespace particle {

  template <typename Ptc >
  void instant_produce_pairs( Real dt, Rng& rng, Ptc& ptc, Ptc& electron, Ptc& positron );

  template <typename Ptc >
  void photon_produce_pairs( Real dt, Ptc& photon, Ptc& electron, Ptc& positron, Rng& rng );

  template <typename Ptc >
  void produce_photons( Ptc& ptc, Ptc& photon, Rng& rng );
}

#endif
