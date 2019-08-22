#ifndef _APT_CSI_HPP_
#define _APT_CSI_HPP_

namespace apt {
  // comma-separated-integer
  template < typename T >
  std::string csi( T x ) {
    T r = x % 1000;
    x /= 1000;
    return x != 0 ? csi(x) + "," + std::to_string(std::abs(r)) : std::to_string(r);
  }
}


#endif
