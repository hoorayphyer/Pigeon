#include "pic.hpp"
#include "simulator/builder_impl.hpp"
namespace pic {
template class SimulationBuilder<DGrid, real_t, particle::Specs, real_j_t,
                                 real_export_t>;
}
