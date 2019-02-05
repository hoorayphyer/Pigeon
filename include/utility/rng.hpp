#ifndef _UTIL_RNG_HPP_
#define _UTIL_RNG_HPP_

#include <random>

// TODO use template ?
namespace util {
  template < typename T>
  class Rng {
  private:
    std::default_random_engine _engine;
    std::uniform_real_distribution<T> _uniform_dist;
    std::normal_distribution<T> _gaussian_dist;

  public:
    inline void set_seed( unsigned long int seed ) {
      _engine.seed( seed );
    }

    // uniform distribution between [0.0, 1.0)
    inline auto uniform() {
      return _uniform_dist(_engine);
    }

    // standard gaussian distribution
    inline auto gaussian() {
      return _gaussian_dist(_engine);
    }

  };
}


#endif
