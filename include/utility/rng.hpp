#ifndef _RNG_HPP_
#define _RNG_HPP_

#include "traits.hpp"
#include <random>

class Rng {
private:
  std::default_random_engine _engine;
  std::uniform_real_distribution<double> _uniform_dist;
  std::normal_distribution<double> _gaussian_dist;

public:
  inline void set_seed( unsigned long int seed ) {
    _engine.seed( seed );
  }

  // uniform distribution between [0.0, 1.0)
  inline double uniform() {
    return _uniform_dist(_engine);
  }

  // standard gaussian distribution
  inline double gaussian() {
    return _gaussian_dist(_engine);
  }

};

#endif // ----- end of _RNG_H_
