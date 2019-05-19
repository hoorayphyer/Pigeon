#include "field/updater_impl.hpp"
#include "pic.hpp"
// only applicable in 2D log spherical
namespace field {
  template struct Updater<pic::real_t, pic::DGrid, pic::real_j_t>;
}
