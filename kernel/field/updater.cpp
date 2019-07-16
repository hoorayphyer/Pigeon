#include "field/updater_impl.hpp"
#include "pic.hpp"

// only applicable in 2D log spherical
namespace field::ofs {
  int magnetic_pole = 2;
  apt::array<int,4> indent {};
  double damping_rate = 10.0;
}

namespace field {
  template struct Updater<pic::real_t, pic::DGrid, pic::real_j_t>;
}
