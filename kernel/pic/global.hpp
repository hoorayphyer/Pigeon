#ifndef _PIC_GLOBAL_HPP_
#define _PIC_GLOBAL_HPP_

#include "particle/forces.hpp"
#include "particle/scattering.hpp"

// TODO use just global function to generate, this will save a lot of template
// contaminations. Also, check ForceGen::operator() for memory leak
namespace particle {
template <typename Real, template <typename> class Specs,
          template <typename, template <typename> class> class Ptc_t>
ForceGen<Real, Specs, Ptc_t> force_gen;

template <typename Real, template <typename> class Specs>
ScatGen<Real, Specs> scat_gen;
}  // namespace particle

#endif
