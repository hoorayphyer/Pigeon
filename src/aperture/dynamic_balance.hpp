#ifndef _DYNAMIC_BALANCE_HPP_
#define _DYNAMIC_BALANCE_HPP_

#include "particle/array.hpp"
#include "particle/map.hpp"

namespace aperture {
  void detailed_balance();

  // NOTE fields are not taken care of during dynamic_adjust, so data such as pair creation rate on each ensemble is simply lost. The solution is to do dynamic_ajust always afeter a data export, which is reset that kind of data anyway.
  void dynamic_adjust();
}

#endif

