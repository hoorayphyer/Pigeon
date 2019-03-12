#ifndef  _APT_PAIR_HPP_
#define  _APT_PAIR_HPP_

namespace apt {
  template < typename T >
  struct pair {
    T _data[2];

    constexpr T operator[] ( bool lr ) const noexcept { return _data[lr];}
    constexpr T& operator[] ( bool lr ) noexcept { return _data[lr];}
  };
}

inline constexpr bool LFT = false;
inline constexpr bool RGT = true;

#endif
