#ifndef _UTIL_RNG_HPP_
#define _UTIL_RNG_HPP_

#include <random>

namespace util {
template <typename T>
class Rng {
 private:
  std::default_random_engine _engine{};
  std::uniform_real_distribution<T> _uniform_dist{};
  std::normal_distribution<T> _gaussian_dist{};

 public:
  inline void set_seed(unsigned long int seed) { _engine.seed(seed); }

  // uniform distribution between [0.0, 1.0)
  inline auto uniform() { return _uniform_dist(_engine); }

  inline auto uniform(T lb, T ub) {
    return lb + _uniform_dist(_engine) * (ub - lb);
  }

  // standard gaussian distribution
  inline auto gaussian() { return _gaussian_dist(_engine); }

  inline auto gaussian(T mu, T sig) {
    return mu + _gaussian_dist(_engine) * sig;
  }
};
}  // namespace util

#endif
