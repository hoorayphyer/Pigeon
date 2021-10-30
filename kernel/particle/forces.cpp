#include "particle/forces_impl.hpp"
#include "pic.hpp"

using namespace pic;
namespace particle {
template struct Force<real_t, Specs>;
}

namespace particle::force {
template void lorentz(vParticle<real_t, Specs>&, real_t,
                      const apt::Vec<real_t, Specs<real_t>::Dim>&,
                      const apt::Vec<real_t, Specs<real_t>::Dim>&, real_t);

template void lorentz_exact(vParticle<real_t, Specs>&, real_t,
                            const apt::Vec<real_t, Specs<real_t>::Dim>&,
                            const apt::Vec<real_t, Specs<real_t>::Dim>&,
                            real_t);
}  // namespace particle::force
