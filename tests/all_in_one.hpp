#ifndef _TEST_ALL_IN_ONE_HPP_
#define _TEST_ALL_IN_ONE_HPP_

#include "apt/print.hpp"
#include "catch2/catch.hpp"
#include <random>
#include <ctime>
#include <cmath>
#include <iostream>
namespace aio {
  using eng_t = std::default_random_engine;
  auto now() { return std::time(0); }

  template < typename T = int >
  using unif_int = std::uniform_int_distribution<T>; // [lower,upper]

  template < typename T = double >
  using unif_real = std::uniform_real_distribution<T>;
}

#endif
