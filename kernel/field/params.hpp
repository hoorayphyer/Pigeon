#ifndef _FIELDUPDATER_PARAMS_HPP_
#define _FIELDUPDATER_PARAMS_HPP_
#include "apt/array.hpp"

namespace field :: ofs {
  extern int magnetic_pole; // 1 for mono-, 2 for di-
  extern apt::array<int,4> indent;
  extern double damping_rate;
}

#endif
