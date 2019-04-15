#include "particle/forces.cpp"
#include "particle/virtual_particle.hpp"
#include "traits.hpp"

using namespace traits;

namespace particle::force {
  using Ptc = vParticle<real_t, DPtc, ptc_state_t>;

  template
  void lorentz<Ptc>( Ptc&, ts::Real<Ptc>, const ts::Vec<Ptc>&, const ts::Vec<Ptc>&, ts::Real<Ptc> );

  template
  void landau0<Ptc>( Ptc&, ts::Real<Ptc>, const ts::Vec<Ptc>&, const ts::Vec<Ptc>&, ts::Real<Ptc> );
}
