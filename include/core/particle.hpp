#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_

#include "types.hpp"

// single particle
class Particle {
private:
  std::array<Real, 7> _internal;
public:
  constexpr Particle( Real q1, Real q2, Real q3, Real p1, Real p2, Real p3, Real state)
    : _internal{q1, q2, q3, p1, p2, p3, state} {}

  constexpr auto q() const {
    return std::forward_as_tuple( _internal[0], _internal[1], _internal[2] );
  }
  constexpr auto q() {
    return std::forward_as_tuple( _internal[0], _internal[1], _internal[2] );
  }

  constexpr auto p() const {
    return std::forward_as_tuple( _internal[3], _internal[4], _internal[5] );
  }
  constexpr auto p() {
    return std::forward_as_tuple( _internal[3], _internal[4], _internal[5] );
  }

  constexpr const auto& state() const {
    return _internal[6];
  }
  constexpr auto& state() {
    return _internal[6];
  }

};

#endif
