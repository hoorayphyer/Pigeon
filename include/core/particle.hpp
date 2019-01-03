#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_

#include "types.hpp"
#include <array>
#include <vector>

template < int Dim_Ptc >
struct Particles {
  static constexpr DPtc = Dim_Ptc;

  std::array<std::vector<Real>, Dim_Ptc> q;
  std::array<std::vector<Real>, Dim_Ptc> p;
  std::vector<encoded_bits_t> state;
};

#endif
