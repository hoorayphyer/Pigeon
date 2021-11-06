#include "pic.hpp"
#include "simulator/simulator_impl.hpp"
namespace pic {
template class Simulator<DGrid, real_t, particle::Specs, real_j_t,
                         real_export_t>;
}
