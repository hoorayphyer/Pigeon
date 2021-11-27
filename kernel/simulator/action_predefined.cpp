#include "pic.hpp"
#include "simulator/action_predefined_impl.hpp"
namespace particle {
template class Migrator<pic::DGrid, pic::real_t, particle::Specs,
                        pic::real_j_t>;
}
