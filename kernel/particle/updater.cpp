#include "particle/updater_impl.hpp"
#include "pic.hpp"
using namespace pic;
namespace particle {
template class Updater<DGrid, real_t, Specs, ShapeF, real_j_t>;
}
