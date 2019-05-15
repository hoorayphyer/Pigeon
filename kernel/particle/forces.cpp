#include "forces_impl.hpp"
#include "pic.hpp"

using namespace pic;

namespace particle::force {
  template
  void lorentz( vParticle<real_t,Specs>&, real_t, const apt::Vec<real_t,Specs<real_t>::Dim>&, const apt::Vec<real_t,Specs<real_t>::Dim>&, real_t );
}
