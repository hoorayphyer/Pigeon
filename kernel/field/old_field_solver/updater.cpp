#include "field/old_field_solver/updater_impl.hpp"
#include "pic.hpp"

namespace field {
  template struct OldSolve<pic::real_t, pic::DGrid, pic::real_j_t>;
}
