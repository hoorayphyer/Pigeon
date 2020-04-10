#include "particle/annihilation_impl.hpp"
#include "pic.hpp"

using namespace pic;

namespace particle {
  template class Annihilator< DGrid, real_t, Specs, ShapeF, real_j_t >;
}
