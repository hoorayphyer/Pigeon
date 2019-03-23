#ifndef  _APT_PAIR_HPP_
#define  _APT_PAIR_HPP_

namespace apt {
  template < typename T >
  struct pair {
    T lft;
    T rgt;

    constexpr T operator[] ( bool lr ) const noexcept { return lr ? rgt : lft;}
    constexpr T& operator[] ( bool lr ) noexcept { return lr ? rgt : lft;}
  };
}

inline constexpr bool LFT = false;
inline constexpr bool RGT = true;

#endif
