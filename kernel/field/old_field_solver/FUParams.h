#ifndef _FUPARAMS_H_
#define _FUPARAMS_H_

#include "Predefs.h"
#include <array>
#include <unordered_map>

struct FBC {
  typedef Scalar (*BCFunc_split_t) (Scalar t);
  typedef Scalar (*BCFunc_split_x) (Scalar x1, Scalar x2, Scalar x3);

  FieldBCType type;
  int indent;

  // for rotating conductor
  BCFunc_split_t ft;
  BCFunc_split_x B1;
  BCFunc_split_x B2;
  BCFunc_split_x B3;
  BCFunc_split_x E1;
  BCFunc_split_x E2;
  BCFunc_split_x E3;

  // for damping
  Scalar damping_rate;
  // Scalar thickness; // physical value
  // std::function<Scalar(Scalar, Scalar)> profile;
};


// TODO listen on locale
struct FUParams {
  std::array<int, 3> neighbor_left{};
  std::array<int, 3> neighbor_right{};
  Grid grid;

  std::array<bool, 6> is_at_boundary() const {
    std::array<bool, 6> res;
    for ( int i = 0; i < 3; ++i ) {
      res[2*i] = ( neighbor_left[i] == NEIGHBOR_NULL );
      res[2*i + 1] = ( neighbor_right[i] == NEIGHBOR_NULL );
    }
    return res;
  }

  std::unordered_map< BoundaryPosition, FBC > fieldBC;
};


#endif
