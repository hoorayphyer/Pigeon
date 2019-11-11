#ifndef _APT_CSI_HPP_
#define _APT_CSI_HPP_
#include <cstdio>

namespace apt {
  // comma-separated-integer
  template < typename T >
  std::string csi( T x ) {
    int r = x % 1000;
    x /= 1000;
    if ( x != 0 ) {
      char str[4];
      std::snprintf(str, 4, "%03d", std::abs(r));
      return csi(x) + "," + str;
    } else {
      return std::to_string(r);
    }
  }
}


#endif
