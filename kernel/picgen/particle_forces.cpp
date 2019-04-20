#include "particle/forces.cpp"
#include "pic.hpp"

using namespace pic;

namespace particle::force {
  template
  void lorentz( vParticle<real_t,Specs>&, real_t, const Vec<real_t,Specs>&, const Vec<real_t,Specs>&, real_t );
}
