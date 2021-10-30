#include "field/log_spherical_solver/updater_impl.hpp"
#include "pic.hpp"

namespace field {
template struct LogSphericalSolver<pic::real_t, pic::DGrid, pic::real_j_t>;
}  // namespace field
